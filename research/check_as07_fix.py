"""Check AS-07 candidate patch against the complete upstream Pyrolite grid."""
from check_as04_as07 import read, section, compile_run, COMMON
from pathlib import Path
import subprocess
import tempfile

source=read("source/material_model/utilities.cc")
with tempfile.TemporaryDirectory() as d:
    root=Path(d)
    p=root/"source/material_model/utilities.cc"
    p.parent.mkdir(parents=True)
    p.write_text(source)
    patch=Path(__file__).with_name("as07.patch").resolve()
    subprocess.run(["git","apply","--check","--unidiff-zero",str(patch)],cwd=root,check=True)
    subprocess.run(["git","apply","--unidiff-zero",str(patch)],cwd=root,check=True)
    fixed=p.read_text()

def methods(s):
    return (section(s,"        double\n        MaterialLookup::value (","        std::array<double,2>\n        MaterialLookup::get_pT_steps")+
            section(s,"        double\n        MaterialLookup::get_nT","        HeFESToReader::HeFESToReader")+
            section(s,"        double\n        MaterialLookup::dRhodp","        unsigned int")+
            section(s,"        double\n        MaterialLookup::dHdT","        std::array<std::pair<double, unsigned int>,2>"))

data=read("data/material-model/steinberger/pyr_MS95_with_volume_fractions_lo_res.dat").splitlines()
pmin=float(data[4])*1e5; dp=float(data[5])*1e5; np=int(data[6])
tmin=float(data[8]); dt=float(data[9]); nt=int(data[10])
assert (pmin,dp,np,tmin,dt,nt)==(1e9,2e9,68,400,100,37)
rho=[[0.]*np for _ in range(nt)]
h=[[0.]*np for _ in range(nt)]
for line in data[13:]:
    if line.strip():
        row=list(map(float,line.split()[:8]))
        i=round((row[1]-tmin)/dt); j=round((row[0]*1e5-pmin)/dp)
        rho[i][j]=row[2]; h[i][j]=row[7]
assert all(x>0 for row in rho for x in row)
def initializer(t):
    return "{"+",".join("{"+",".join(map(str,row))+"}" for row in t)+"}"

pre=r"""
template<int dim,class T> struct Table {
 std::vector<std::vector<T>> a;
 auto n_rows()const{return a.size();}
 auto n_cols()const{return a[0].size();}
 const auto& operator[](unsigned i)const{return a.at(i);}
};
struct MaterialLookup {
 double min_temp=400,max_temp=4000,delta_temp=100,min_press=1e9,max_press=135e9,delta_press=2e9;
 unsigned n_temperature=37,n_pressure=68;
 bool interpolation=true;
 Table<2,double> density_values,enthalpy_values;
 double get_nT(double)const;
 double get_np(double)const;
 double value(double,double,const Table<2,double>&,bool)const;
 unsigned int value(double,double,const Table<2,unsigned>&)const;
 double dRhodp(double,double)const;
 double dHdT(double,double)const;
 double dHdp(double,double)const;
};
"""
post=r"""
bool near(double a,double b){return std::abs(a-b)<=1e-10*std::max(std::abs(b),1e-12);}
int main(){
 MaterialLookup m;
 m.density_values.a=RHO_DATA;
 m.enthalpy_values.a=H_DATA;
 auto &a=m.density_values.a;
 auto &h=m.enthalpy_values.a;
 Table<2,unsigned> categorical;
 categorical.a.resize(37,std::vector<unsigned>(68));
 unsigned knots=0,midpoints=0,derivatives=0;
 for(unsigned i=0;i<37;++i)for(unsigned j=0;j<68;++j){
  double t=400+100*i,p=1e9+2e9*j;
  categorical.a[i][j]=i*68+j;
  assert(near(m.value(t,p,m.density_values,true),a[i][j]));
  assert(near(m.value(t,p,m.density_values,false),a[i][j]));
  assert(m.value(t,p,categorical)==i*68+j);
  ++knots;
  if(i<36 && j<67){
   double expected=(a[i][j]+a[i+1][j]+a[i][j+1]+a[i+1][j+1])/4;
   assert(near(m.value(t+50,p+1e9,m.density_values,true),expected));
   ++midpoints;
  }
  unsigned k=std::min(j,66u),l=std::min(i,35u);
  assert(near(m.dRhodp(t,p),(a[i][k+1]-a[i][k])/2e9));
  assert(near(m.dHdp(t,p),(h[i][k+1]-h[i][k])/2e9));
  assert(near(m.dHdT(t,p),(h[l+1][j]-h[l][j])/100));
  derivatives+=3;
 }
 for(bool interpol:{false,true}){
  m.interpolation=interpol;
  assert(near(m.dRhodp(2000,134e9),(a[16][67]-a[16][66])/2e9));
  assert(near(m.dRhodp(2000,135e9),(a[16][67]-a[16][66])/2e9));
  assert(near(m.dHdT(3950,51e9),(h[36][25]-h[35][25])/100));
  assert(near(m.dHdp(2000,134e9),(h[16][67]-h[16][66])/2e9));
  assert(near(m.value(5000,150e9,m.density_values,interpol),a.back().back()));
  assert(near(m.value(100,0,m.density_values,interpol),a[0][0]));
 }
 std::cout<<std::setprecision(12)<<"134 GPa density="<<m.value(2000,134e9,m.density_values,true)<<" derivative="<<m.dRhodp(2000,134e9)<<"\n";
 std::cout<<"135 GPa density="<<m.value(2000,135e9,m.density_values,true)<<"\n";
 std::cout<<"3950 K density="<<m.value(3950,51e9,m.density_values,true)<<"\n";
 std::cout<<"AS-07 FIX PASS: "<<knots<<" grid knots (3 lookup modes), "<<midpoints<<" cell midpoints, "<<derivatives<<" finite differences, 12 edge controls\n";
}
"""
compile_run("AS-07 candidate",COMMON+pre+methods(fixed)+post.replace("RHO_DATA",initializer(rho)).replace("H_DATA",initializer(h)))
