#ifndef Utilities_bin_h
#define Utilities_bin_h

#include <iostream>
#include <tuple>
#include <string>
#include <set>

// a simple template class to store a bin and overload the equality operator

// define a few common uses of the template class
template <typename T> class bin : public std::tuple<std::string,T,T,T,T>
{
 public:
  bin() : std::tuple<std::string,T,T,T,T>() {};
  bin(std::string n, T a, T b, T c=0, T d=0) : std::tuple<std::string,T,T,T,T>(n,a,b,c,d) {};
  std::string name() const { return std::get<0>(*this); }
  T low()   const { return std::get<1>(*this); }
  T high()  const { return std::get<2>(*this); }
  T mean()  const { return std::get<3>(*this); }
  T width() const { return std::get<4>(*this); }
  void print() const {
    std::string l = (this->name()+" = [ "+std::to_string(this->mean())+" , "+std::to_string(this->width())+" , "+std::to_string(this->low())+" , "+std::to_string(this->high())+" ]");
    std::cout << l << std::endl;
  }
  bool operator < (const bin& ref) const 
  {
    const std::tuple<std::string,T,T> a(this->name(), this->low(), this->high());
    const std::tuple<std::string,T,T> b(ref.name(), ref.low(), ref.high());
    return (bool)(a < b);
  }
};
typedef bin<float>  binF;

// associate a set of bins to make an analysis bin
class anabin : public std::set<binF>
{
 public:
  anabin() : std::set<binF>() {};
  anabin(const std::set<binF>& set) : std::set<binF>(set) {};
  binF getbin (const std::string& name) const {
    for (const auto& s : *this) { if (s.name()==name) { return s; } };
    return binF("", -99., -99., -99., -99.);
  };
  void setbin (const binF& bin) { this->insert(bin); };
  void setbin (const std::string& n, const float& a, const float& b, const float& c=-99., const float& d=-99.) {
    const auto& mean  = (c!=-99. ? c : (b+a)/2.);
    const auto& width = (d!=-99. ? d : (b-a)/2.);
    this->insert(binF(n,a,b,mean,width));
  };
  void print() const {
    std::string l = "";
    for (const auto& s : *this) { l += (s.name()+" = [ "+std::to_string(s.mean())+" , "+std::to_string(s.width())+" , "+std::to_string(s.low())+" , "+std::to_string(s.high())+" ] , "); }
    if (l.rfind("] , ")!=std::string::npos) { l.erase(l.rfind(" , "), 3); }
    std::cout << l << std::endl;
  }
};

#endif // #ifndef bin_h
