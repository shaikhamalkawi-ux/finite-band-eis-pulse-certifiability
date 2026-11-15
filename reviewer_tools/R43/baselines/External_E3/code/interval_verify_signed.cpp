// Independent directed-rounding MPFR verifier for the positive-RC continuum dual.
// No optimizer or mesh-convergence assertion enters acceptance.
// Build: g++ -O2 -std=c++17 interval_verify.cpp -Wl,-l:libmpfr.so.6 -lgmp -o interval_verify
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <gmp.h>
#if __has_include(<mpfr.h>)
#include <mpfr.h>
#else
// Minimal public MPFR 4.x ABI declarations for a runtime without development headers.
// On normal systems use the official mpfr.h via the branch above.
extern "C" {
typedef long mpfr_prec_t; typedef long mpfr_exp_t;
typedef struct {mpfr_prec_t _mpfr_prec;int _mpfr_sign;mpfr_exp_t _mpfr_exp;mp_limb_t *_mpfr_d;} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1]; typedef __mpfr_struct *mpfr_ptr; typedef const __mpfr_struct *mpfr_srcptr;
typedef enum {MPFR_RNDN=0,MPFR_RNDZ,MPFR_RNDU,MPFR_RNDD,MPFR_RNDA,MPFR_RNDF} mpfr_rnd_t;
void mpfr_init2(mpfr_ptr,mpfr_prec_t);void mpfr_clear(mpfr_ptr);
int mpfr_set_d(mpfr_ptr,double,mpfr_rnd_t);int mpfr_set(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_add(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_sub(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_mul(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_div(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_neg(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_abs(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_exp(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_expm1(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_const_pi(mpfr_ptr,mpfr_rnd_t);
int mpfr_cmp(mpfr_srcptr,mpfr_srcptr);int mpfr_cmp_si(mpfr_srcptr,long);int mpfr_number_p(mpfr_srcptr);
double mpfr_get_d(mpfr_srcptr,mpfr_rnd_t);const char* mpfr_get_version(void);
}
#endif
static int PREC=160;
struct Iv {
 mpfr_t lo,hi;
 Iv(double v=0){mpfr_init2(lo,PREC);mpfr_init2(hi,PREC);mpfr_set_d(lo,v,MPFR_RNDD);mpfr_set_d(hi,v,MPFR_RNDU);}
 Iv(double a,double b):Iv(){if(a>b)throw std::runtime_error("bad interval");mpfr_set_d(lo,a,MPFR_RNDD);mpfr_set_d(hi,b,MPFR_RNDU);}
 Iv(const Iv& o):Iv(){mpfr_set(lo,o.lo,MPFR_RNDD);mpfr_set(hi,o.hi,MPFR_RNDU);}
 Iv& operator=(const Iv&o){if(this!=&o){mpfr_set(lo,o.lo,MPFR_RNDD);mpfr_set(hi,o.hi,MPFR_RNDU);}return *this;}
 ~Iv(){mpfr_clear(lo);mpfr_clear(hi);}
 bool positive()const{return mpfr_number_p(lo)&&mpfr_number_p(hi)&&mpfr_cmp_si(lo,0)>0;}
 bool nonnegative()const{return mpfr_number_p(lo)&&mpfr_number_p(hi)&&mpfr_cmp_si(lo,0)>=0;}
 double lower()const{return mpfr_get_d(lo,MPFR_RNDD);} double upper()const{return mpfr_get_d(hi,MPFR_RNDU);}
};
Iv operator+(const Iv&a,const Iv&b){Iv c;mpfr_add(c.lo,a.lo,b.lo,MPFR_RNDD);mpfr_add(c.hi,a.hi,b.hi,MPFR_RNDU);return c;}
Iv operator-(const Iv&a,const Iv&b){Iv c;mpfr_sub(c.lo,a.lo,b.hi,MPFR_RNDD);mpfr_sub(c.hi,a.hi,b.lo,MPFR_RNDU);return c;}
Iv operator-(const Iv&a){Iv c;mpfr_neg(c.lo,a.hi,MPFR_RNDD);mpfr_neg(c.hi,a.lo,MPFR_RNDU);return c;}
Iv operator*(const Iv&a,const Iv&b){
 Iv c,t;
 mpfr_mul(c.lo,a.lo,b.lo,MPFR_RNDD);mpfr_mul(c.hi,a.lo,b.lo,MPFR_RNDU);
 mpfr_srcptr as[2]={a.lo,a.hi},bs[2]={b.lo,b.hi};
 for(int i=0;i<2;i++)for(int j=0;j<2;j++){
  mpfr_mul(t.lo,as[i],bs[j],MPFR_RNDD);if(mpfr_cmp(t.lo,c.lo)<0)mpfr_set(c.lo,t.lo,MPFR_RNDD);
  mpfr_mul(t.hi,as[i],bs[j],MPFR_RNDU);if(mpfr_cmp(t.hi,c.hi)>0)mpfr_set(c.hi,t.hi,MPFR_RNDU);
 }return c;
}
Iv inv(const Iv&a){if(mpfr_cmp_si(a.lo,0)<=0&&mpfr_cmp_si(a.hi,0)>=0)throw std::runtime_error("division through zero");Iv c,one(1);mpfr_div(c.lo,one.lo,a.hi,MPFR_RNDD);mpfr_div(c.hi,one.hi,a.lo,MPFR_RNDU);return c;}
Iv operator/(const Iv&a,const Iv&b){return a*inv(b);}
Iv expi(const Iv&a){Iv c;mpfr_exp(c.lo,a.lo,MPFR_RNDD);mpfr_exp(c.hi,a.hi,MPFR_RNDU);return c;}
Iv expm1i(const Iv&a){Iv c;mpfr_expm1(c.lo,a.lo,MPFR_RNDD);mpfr_expm1(c.hi,a.hi,MPFR_RNDU);return c;}
Iv absi(const Iv&a){if(mpfr_cmp_si(a.lo,0)>=0)return a;if(mpfr_cmp_si(a.hi,0)<=0)return -a;Iv c,t;mpfr_abs(c.hi,a.lo,MPFR_RNDU);mpfr_abs(t.hi,a.hi,MPFR_RNDU);if(mpfr_cmp(t.hi,c.hi)>0)mpfr_set(c.hi,t.hi,MPFR_RNDU);return c;}
Iv sq(const Iv&a){if(mpfr_cmp_si(a.lo,0)>=0){Iv c;mpfr_mul(c.lo,a.lo,a.lo,MPFR_RNDD);mpfr_mul(c.hi,a.hi,a.hi,MPFR_RNDU);return c;}if(mpfr_cmp_si(a.hi,0)<=0)return sq(-a);Iv b=absi(a);return sq(b);}
Iv pi(){Iv v;mpfr_const_pi(v.lo,MPFR_RNDD);mpfr_const_pi(v.hi,MPFR_RNDU);return v;}
struct Model {
 int n,p;double iT,eps,meas,alpha,T;std::vector<Iv>w,yr,yi,zr,zi,aa,bb,u,delay,jump;Iv Q,M1,last;
 Model(std::string path){std::ifstream in(path);if(!in)throw std::runtime_error("input missing");in>>n>>p>>T>>iT>>eps>>meas>>alpha;
  if(!in || n<1 || n>10000 || p<1 || p>100000 || !std::isfinite(T) || T<=0 || !std::isfinite(iT) || !std::isfinite(eps) || eps<=0 || !std::isfinite(meas) || !std::isfinite(alpha) || alpha<=0) throw std::runtime_error("invalid header");
  double previous_end=0;
  for(int k=0;k<n;k++){double f,r,i,yrr,yii;in>>f>>r>>i>>yrr>>yii;if(!in || !std::isfinite(f) || f<=0 || !std::isfinite(r) || !std::isfinite(i) || !std::isfinite(yrr) || !std::isfinite(yii))throw std::runtime_error("invalid frequency row");w.push_back(Iv(2)*pi()*Iv(f));zr.emplace_back(r);zi.emplace_back(i);yr.emplace_back(yrr);yi.emplace_back(yii);}
  for(int j=0;j<p;j++){double a,b,v;in>>a>>b>>v;if(!in || !std::isfinite(a) || !std::isfinite(b) || !std::isfinite(v) || !(a<b&&a==previous_end&&b<=T))throw std::runtime_error("bad waveform");previous_end=b;aa.emplace_back(a);bb.emplace_back(b);u.emplace_back(v);Q=Q+Iv(v)*(Iv(b)-Iv(a));M1=M1+absi(Iv(v))*(sq(Iv(T)-Iv(a))-sq(Iv(T)-Iv(b)))/Iv(2);}
  if(!in || previous_end!=T)throw std::runtime_error("incomplete waveform");std::string extra;if(in>>extra)throw std::runtime_error("trailing input");last=u.back();
  delay.push_back(Iv(T)-aa[0]);jump.push_back(-u[0]);
  for(int j=1;j<p;j++){Iv d=u[j-1]-u[j];if(d.lower()!=0||d.upper()!=0){delay.push_back(Iv(T)-aa[j]);jump.push_back(d);}}
 }
 Iv residual(const Iv&x,int deriv=0)const{
  Iv tau=expi(x),total=deriv==0?last:Iv(0);
  for(int k=0;k<n;k++){
   Iv v=w[k]*tau,z=sq(v),den=Iv(1)+z;
   if(deriv==0)total=total+(yr[k]-yi[k]*v)/den;
   if(deriv==1)total=total+yr[k]*Iv(-2)*z/sq(den)+yi[k]*v*(z-Iv(1))/sq(den);
   if(deriv==2)total=total+(yr[k]*Iv(4)*z*(z-Iv(1))+yi[k]*v*(-sq(z)+Iv(6)*z-Iv(1)))/(sq(den)*den);
  }
  for(size_t j=0;j<jump.size();j++){
   Iv q=delay[j]/tau,e=expi(-q);
   if(deriv==0)total=total+jump[j]*e;
   if(deriv==1)total=total+jump[j]*q*e;
   if(deriv==2)total=total+jump[j]*q*(q-Iv(1))*e;
  }
  return total;
 }
 Iv objective()const{Iv c;for(int k=0;k<n;k++)c=c+yr[k]*zr[k]+yi[k]*zi[k]+Iv(alpha)*Iv(eps)*(absi(yr[k])+absi(yi[k]));return -c;}
 std::vector<Iv> finite()const{Iv r(iT),l,c(Q);for(int k=0;k<n;k++){r=r+yr[k];l=l+yi[k]*w[k];c=c-yi[k]/w[k];}return {r,l,c};}
 Iv smalltail(double xl)const{
  Iv tau=expi(Iv(xl)),s(last);
  for(int k=0;k<n;k++){s=s+yr[k]-absi(yr[k])*sq(w[k]*tau)-absi(yi[k])*w[k]*tau;}
  for(size_t j=0;j<jump.size();j++)s=s-absi(jump[j])*expi(-delay[j]/tau);
  return s;
 }
 Iv largetail(double xh)const{
  // r(1/z)/z = H(z). For signed u, |exp(-q)-1|<=q; M1 integrates |u|*(T-s).
  Iv z=expi(Iv(-xh)),s=finite()[2]-z*M1;
  for(int k=0;k<n;k++)s=s-absi(yr[k])*z/sq(w[k])-absi(yi[k])*sq(z)/(sq(w[k])*w[k]);
  return s;
 }
};
struct Bin {double a,b;int depth;};
int main(int argc,char**argv){try{
 if(argc<3){std::cerr<<"usage: verifier input.txt leaves.csv [bits=160] [mode=2] [shift=0]\n";return 2;}
 if(argc>3)PREC=std::stoi(argv[3]);int mode=argc>4?std::stoi(argv[4]):2;int shift=argc>5?std::stoi(argv[5]):0;
 if(PREC<64 || PREC>4096 || (mode!=1&&mode!=2) || shift<0 || shift>32)throw std::runtime_error("invalid verifier options");Model m(argv[1]);auto fs=m.finite();Iv bound=m.objective(),gap=bound-Iv(m.meas);
 double xl=-24-shift,xh=8+shift;Iv small=m.smalltail(xl),large=m.largetail(xh);
 bool pre=fs[0].nonnegative()&&fs[1].nonnegative()&&fs[2].nonnegative()&&small.positive()&&large.positive()&&gap.positive();
 std::ofstream out(argv[2]);out<<"x_lo,x_hi,residual_lower,method,depth\n"<<std::setprecision(17);
 std::vector<Bin> stack;for(double a=xl;a<xh;a+=1)stack.push_back({a,std::min(xh,a+1),0});
 long evals=0,accepted=0,nn=0,mv=0,tt=0;double minlb=1e300;int maxdepth=0;bool compact=true;
 if(pre)while(!stack.empty()){
  auto b=stack.back();stack.pop_back();Iv X(b.a,b.b);evals++;Iv r=m.residual(X);int method=0;
  if(!r.positive()){
   double mid=(b.a+b.b)/2,h=(b.b-b.a)/2;Iv center=m.residual(Iv(mid));
   if(mode==1){r=center+Iv(-h,h)*m.residual(X,1);method=1;}
   else {Iv second=m.residual(X,2);r=center+Iv(-h,h)*m.residual(Iv(mid),1)+Iv(0,h*h/2)*second;method=2;}
  }
  if(r.positive()){accepted++;minlb=std::min(minlb,r.lower());maxdepth=std::max(maxdepth,b.depth);if(method==0)nn++;else if(method==1)mv++;else tt++;out<<b.a<<","<<b.b<<","<<r.lower()<<","<<method<<","<<b.depth<<"\n";}
  else {if(b.depth>=28||evals>200000){compact=false;out<<b.a<<","<<b.b<<","<<r.lower()<<",FAIL,"<<b.depth<<"\n";break;}double mid=(b.a+b.b)/2;stack.push_back({mid,b.b,b.depth+1});stack.push_back({b.a,mid,b.depth+1});}
 } else compact=false;
 bool pass=pre&&compact&&stack.empty();
 std::cout<<std::setprecision(17)<<"{\"passed\":"<<(pass?"true":"false")<<",\"precision_bits\":"<<PREC<<",\"mpfr_version\":\""<<mpfr_get_version()<<"\",\"mode\":"<<mode<<",\"xl\":"<<xl<<",\"xh\":"<<xh<<",\"certified_lower_mOhm\":"<<bound.lower()<<",\"certified_lower_upper_round_mOhm\":"<<bound.upper()<<",\"gap_lower_mOhm\":"<<gap.lower()<<",\"finite_R0_lower\":"<<fs[0].lower()<<",\"finite_L_lower\":"<<fs[1].lower()<<",\"finite_c_lower\":"<<fs[2].lower()<<",\"fast_tail_lower\":"<<small.lower()<<",\"slow_tail_H_lower\":"<<large.lower()<<",\"compact_residual_lower\":"<<(accepted?minlb:0)<<",\"accepted_leaves\":"<<accepted<<",\"evaluated_bins\":"<<evals<<",\"max_depth\":"<<maxdepth<<",\"natural_leaves\":"<<nn<<",\"mean_value_leaves\":"<<mv<<",\"taylor_leaves\":"<<tt<<"}\n";
 return pass?0:1;
 }catch(std::exception&e){std::cerr<<e.what()<<"\n";return 3;}}
