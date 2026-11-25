"""Compile selected production arithmetic against independent invariants.
Source and table are pinned upstream files. This is not an ASPECT/deal.II integration test.
"""
import pathlib, subprocess, tempfile, urllib.request, re, sys

SHA = "540a201c460759e160968bfc03323db8e7b76367"
ROOT = f"https://raw.githubusercontent.com/geodynamics/aspect/{SHA}/"

def read(path):
    with urllib.request.urlopen(ROOT + path, timeout=40) as r:
        return r.read().decode()

def section(s, begin, end):
    a = s.index(begin)
    return s[a:s.index(end, a)]

def compile_run(name, source):
    with tempfile.TemporaryDirectory() as d:
        p = pathlib.Path(d)
        (p / "test.cc").write_text(source)
        subprocess.run(["c++", "-std=c++17", "-O2", "-Wall", str(p / "test.cc"), "-o", str(p / "test")], check=True)
        print(name, flush=True)
        subprocess.run([str(p / "test")], check=True)

COMMON = r"""
#include <algorithm>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <utility>
#define Assert(c,m) assert(c)
#define AssertThrow(c,m) assert(c)
"""

def as04():
    src = read("source/material_model/reaction_model/katz2003_mantle_melting.cc")
    box = read("source/geometry_model/box.cc")
    old = section(src, "// remove melt that gets close to the surface", "// because depletion is a volume-based")
    fix = old.replace("in.position[i](1) - (this->get_geometry_model().maximal_depth() - extraction_depth)", "extraction_depth - this->get_geometry_model().depth(in.position[i])")
    assert fix != old
    depth = re.search(r"const double d = extents\[dim-1\].*?return std::clamp\(d, 0\., maximal_depth\(\)\);",box,re.S).group()
    cpp = r"""
template<int dim> struct Point {
 std::array<double,dim> x;
 double operator()(unsigned i) const {return x[i];}
};
template<int dim> struct Box {
 std::array<double,dim> extents,box_origin;
 double maximal_depth() const {return extents[dim-1];}
 double depth(Point<dim> position) const {
 const double topo=0;
 DEPTH_BODY
 }
};
template<int dim> struct Model {
 Box<dim> geometry;
 const Box<dim>& get_geometry_model() const {return geometry;}
 double old(Point<dim> point,double extraction_depth,double old_porosity,double porosity_change) const {
 struct {std::vector<Point<dim>> position;} in{{point}};
 unsigned i=0;
 OLD_BODY
 return porosity_change;
 }
 double fixed(Point<dim> point,double extraction_depth,double old_porosity,double porosity_change) const {
 struct {std::vector<Point<dim>> position;} in{{point}};
 unsigned i=0;
 FIX_BODY
 return porosity_change;
 }
};
template<int dim> void test(const char* label, Model<dim> model,Point<dim> p,double expected,bool same,double extraction_depth=1000,double prior=0) {
 double a=model.old(p,extraction_depth,.05,prior);
 double b=model.fixed(p,extraction_depth,.05,prior);
 std::cout<<label<<" depth="<<model.geometry.depth(p)<<" expected="<<expected<<" old="<<a<<" fixed="<<b<<"\n";
 assert(std::abs(b-expected)<1e-12);
 if(same) assert(std::abs(a-expected)<1e-12);
 else assert(std::abs(a-expected)>1e-8);
}
int main() {
 std::cout<<std::setprecision(10);
 Model<2> m2{{{100000,100000},{0,0}}};
 Model<3> m3{{{100000,100000,100000},{0,0,0}}};
 Model<2> shifted{{{100000,100000},{0,-10000}}};
 test<2>("2D control",m2,{{50000,99500}},-.025,true);
 test<3>("3D y=10km",m3,{{50000,10000,99500}},-.025,false);
 test<3>("3D y=50km",m3,{{50000,50000,99500}},-.025,false);
 test<3>("3D y=99.5km",m3,{{50000,99500,99500}},-.025,true);
 test<2>("2D translated",shifted,{{50000,89500}},-.025,false);
 test<3>("outside layer",m3,{{50000,50000,98000}},.01,true,1000,.01);
 test<3>("extraction disabled",m3,{{50000,50000,100000}},.01,true,0,.01);
 std::cout<<"AS-04 PASS: 7 cases; 3 counterexamples\n";
}
"""
    compile_run("AS-04", COMMON+cpp.replace("DEPTH_BODY",depth).replace("OLD_BODY",old).replace("FIX_BODY",fix))

def as05():
    s=read("source/material_model/utilities.cc")
    old=section(s,"      void\n      reaction_progress_modify_values","      double phase_average_value")
    rest=section(s,"      double phase_average_value","      template <int dim>\n      PhaseFunctionInputs")
    fixed=old.replace("previous_phase_index);","previous_phase_index - composition_index);").replace("reaction_progress_modify_values (","fixed_modify_values (")
    assert "previous_phase_index - composition_index" in fixed
    pre=r"""
namespace PhaseUtilities { enum PhaseAveragingOperation {arithmetic,logarithmic}; }
std::pair<bool,double> get_reaction_progress_for_phase_transition(const std::vector<double>&,const std::vector<unsigned>&,unsigned);
"""
    post=r"""
void run(const char* label, std::vector<unsigned> n, unsigned map, unsigned comp, std::vector<double> rho, double progress, double expected, bool same) {
 std::vector<double> a=rho,b=rho,phases;
 for(auto v:n)for(unsigned j=0;j<v;j++)phases.push_back(1.);
 for(unsigned c=0;c<n.size();c++){
  reaction_progress_modify_values({progress},{map},n,a,c,PhaseUtilities::arithmetic);
  fixed_modify_values({progress},{map},n,b,c,PhaseUtilities::arithmetic);
 }
 double x=phase_average_value(phases,n,a,comp,PhaseUtilities::arithmetic);
 double y=phase_average_value(phases,n,b,comp,PhaseUtilities::arithmetic);
 std::cout<<label<<" expected="<<expected<<" old="<<x<<" fixed="<<y<<"\n";
 assert(std::abs(y-expected)<1e-9);
 if(same)assert(std::abs(x-expected)<1e-9);else assert(std::abs(x-expected)>1e-9);
}
int main(){
 run("background p=0",{1},0,0,{3000,3540},0,3000,true);
 run("background p=.5",{1},0,0,{3000,3540},.5,3270,true);
 run("crust p=0",{0,1},0,1,{3300,3000,3540},0,3000,false);
 run("crust p=.5",{0,1},0,1,{3300,3000,3540},.5,3270,false);
 run("crust p=1",{0,1},0,1,{3300,3000,3540},1,3540,true);
 run("crust after background transition",{1,1},1,1,{3300,3394,3000,3540},.5,3270,false);
 run("unrelated background",{1,1},1,0,{3300,3394,3000,3540},.5,3394,true);
 run("mapping second transition should not modify first chemical",{0,1,1},1,1,{3300,3000,3540,3235,3372},0,3540,false);
 run("mapping second transition should modify second chemical",{0,1,1},1,2,{3300,3000,3540,3235,3372},0,3235,false);
 std::cout<<"AS-05 PASS: 9 cases; 5 counterexamples\n";
}
"""
    compile_run("AS-05", COMMON+pre+old+fixed+rest+post)

def as06():
    src=read("source/material_model/reactive_fluid_transport.cc")
    body=section(src,"const double total_mass_fraction_water =","// Apply the reactions")
    a=src.index("melt_fractions[q] = in.composition[q][bound_fluid_idx] + mass_fraction_porosity;")
    eq=src[a:src.index(";",a)+1]
    fixed=eq.replace("in.composition[q][bound_fluid_idx]","in.composition[q][bound_fluid_idx] * (1 - mass_fraction_porosity)")
    assert eq != fixed
    cpp=r"""
double compute_bulk_density(double p,double s,double f){return (1-p)*s+p*f;}
double compute_mass_fraction(double v,double d,double b){return v*d/b;}
void test(double volume_fraction_porosity,double water,double solid_density,double fluid_density,double reaction_fraction,bool fix){
 unsigned q=0,bound_fluid_idx=0;
 struct {std::vector<std::vector<double>> composition;} in{{{water}}};
 const double bulk_density=compute_bulk_density(volume_fraction_porosity,solid_density,fluid_density);
 const double mass_fraction_porosity=compute_mass_fraction(volume_fraction_porosity,fluid_density,bulk_density);
 std::vector<double> melt_fractions(1);
 if(fix){ FIX_EQ }else{ OLD_EQ }
 auto eq_free_fluid_fractions=melt_fractions;
 BODY
 double new_p=volume_fraction_porosity+porosity_change;
 double new_w=water+bound_fluid_change;
 double total_after=(new_w*(1-new_p)*solid_density+new_p*fluid_density)/compute_bulk_density(new_p,solid_density,fluid_density);
 double relative=(total_after/total_mass_fraction_water-1)*100;
 std::cout<<"phi="<<volume_fraction_porosity<<" water="<<water<<" fraction="<<reaction_fraction<<" fixed="<<fix<<" total_before="<<total_mass_fraction_water<<" after="<<total_after<<" relative_change_pct="<<relative<<" new_phi="<<new_p<<" new_bound="<<new_w<<"\n";
 if(fix || volume_fraction_porosity==0 || reaction_fraction==0)assert(std::abs(relative)<1e-10);
 else if(reaction_fraction==1)assert(relative>0);
 else assert(std::abs(relative)<1e-10);
}
int main(){
 std::cout<<std::setprecision(12);
 for(double p:{0.,.01,.1,.3}) for(bool f:{false,true})test(p,.01,3000,2000,1.,f);
 for(bool f:{false,true}){test(.1,.01,3000,2000,0.,f);test(.1,.01,3000,2000,.5,f);}
 std::cout<<"AS-06 PASS: 12 executions; 3 conservation counterexamples\n";
}
"""
    compile_run("AS-06", COMMON+cpp.replace("FIX_EQ",fixed).replace("OLD_EQ",eq).replace("BODY",body))

def as07():
    s=read("source/material_model/utilities.cc")
    value=section(s,"        double\n        MaterialLookup::value (","        unsigned int\n        MaterialLookup::value")
    indices=section(s,"        double\n        MaterialLookup::get_nT","        HeFESToReader::HeFESToReader")
    derivative=section(s,"        double\n        MaterialLookup::dRhodp","        unsigned int")
    data=read("data/material-model/steinberger/pyr_MS95_with_volume_fractions_lo_res.dat").splitlines()
    pmin=float(data[4])*1e5; dp=float(data[5])*1e5; np=int(data[6])
    tmin=float(data[8]); dt=float(data[9]); nt=int(data[10])
    rows=[list(map(float,x.split()[:3])) for x in data[13:] if x.strip()]
    assert len(rows)==np*nt
    table=[[0.]*np for i in range(nt)]
    for p,t,rho in rows: table[round((t-tmin)/dt)][round((p*1e5-pmin)/dp)]=rho
    assert all(v>0 for row in table for v in row)
    pre=r"""
template<int dim,class T> struct Table {
 std::vector<std::vector<T>> a;
 auto n_rows()const{return a.size();}
 auto n_cols()const{return a[0].size();}
 const auto& operator[](unsigned i)const{return a.at(i);}
};
struct MaterialLookup {
 double min_temp=400, max_temp=4000,delta_temp=100,min_press=1e9,max_press=135e9,delta_press=2e9;
 bool interpolation=true;
 Table<2,double> density_values;
 double get_nT(double)const;
 double get_np(double)const;
 double value(double,double,const Table<2,double>&,bool)const;
 double dRhodp(double,double)const;
};
"""
    post=r"""
int main(){
 MaterialLookup m;
 m.density_values.a=DATA;
 const double t=2000;
 const unsigned i=16, j=66;
 auto &a=m.density_values.a;
 const double expected=(a[i][j]+a[i][j+1])/2;
 const double old=m.value(t,134e9,m.density_values,true);
 const double last=m.value(t,135e9,m.density_values,true);
 std::cout<<std::setprecision(12)<<"pressure midpoint 134 GPa expected="<<expected<<" actual="<<old<<"; last knot expected="<<a[i][j+1]<<" actual="<<last<<"\n";
 assert(std::abs(old-a[i][j])<1e-9 && std::abs(old-expected)>1e-6);
 assert(std::abs(last-a[i][j])<1e-9 && std::abs(last-a[i][j+1])>1e-6);
 const double expected_dr=(a[i][j+1]-a[i][j])/2e9;
 std::cout<<"last interval dRhodp reference="<<expected_dr<<" actual="<<m.dRhodp(t,134e9)<<"\n";
 assert(expected_dr>0 && m.dRhodp(t,134e9)==0);
 double tmid=m.value(3950,51e9,m.density_values,true);
 double tref=(a[35][25]+a[36][25])/2;
 std::cout<<"temperature midpoint 3950 K expected="<<tref<<" actual="<<tmid<<"\n";
 assert(std::abs(tmid-a[35][25])<1e-9 && std::abs(tmid-tref)>1e-6);
 std::cout<<"interior knot expected="<<a[16][25]<<" actual="<<m.value(2000,51e9,m.density_values,true)<<"\n";
 assert(m.value(2000,51e9,m.density_values,true)==a[16][25]);
 std::cout<<"AS-07 PASS: 5 checks; 4 counterexamples, baseline only\n";
}
"""
    assert (pmin,dp,np,tmin,dt,nt)==(1e9,2e9,68,400,100,37)
    compile_run("AS-07", COMMON+pre+value+indices+derivative+post.replace("DATA","{"+",".join("{"+",".join(map(str,row))+"}" for row in table)+"}"))

if __name__=="__main__":
    tests={"AS-04":as04,"AS-05":as05,"AS-06":as06,"AS-07":as07}
    for name in sys.argv[1:] or tests:
        tests[name]()
