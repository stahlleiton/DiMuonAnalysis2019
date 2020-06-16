#ifndef BIN_H_
#define BIN_H_

#include <iostream>
#include <tuple>
#include <string>
#include <set>
#include <map>

// a simple template class to store a bin and overload the equality operator

// define a few common uses of the template class
template <typename T>
class Bin_t : public std::tuple<std::string,T,T,T,T>
{
 public:
  Bin_t() : std::tuple<std::string,T,T,T,T>() {};
  Bin_t(std::string n, T a, T b, T c=-99, T d=-99) : std::tuple<std::string,T,T,T,T>(n,a,b,c,d) {
    if (c==-99) { std::get<3>(*this) = (b+a)/2.; }
    if (d==-99) { std::get<4>(*this) = (b-a)/2.; }
  };
  const std::string& name() const { return std::get<0>(*this); }
  T low()   const { return std::get<1>(*this); }
  T high()  const { return std::get<2>(*this); }
  T mean()  const { return std::get<3>(*this); }
  T width() const { return std::get<4>(*this); }
  bool valid() const { return this->name()!=""; }
  void print() const {
    std::string l = (this->name()+" = [ "+std::to_string(this->mean())+" , "+std::to_string(this->width())+" , "+std::to_string(this->low())+" , "+std::to_string(this->high())+" ]");
    std::cout << l << std::endl;
  }
  bool operator < (const Bin_t& ref) const 
  {
    const std::tuple<std::string,T,T> a(this->name(), this->low(), this->high());
    const std::tuple<std::string,T,T> b(ref.name(), ref.low(), ref.high());
    return (bool)(a < b);
  }
};
typedef Bin_t<float>  BinF_t;

// associate a set of bins to make an analysis bin
class AnaBin_t
{
 public:
  typedef std::map<std::string, BinF_t> ANABINTYPE;
  class ANABINITER : public ANABINTYPE::iterator
  {
   public:
    ANABINITER() : ANABINTYPE::iterator() {};
    ANABINITER(ANABINTYPE::iterator s) : ANABINTYPE::iterator(s) {};
    BinF_t* operator->() { return (BinF_t* const)&(ANABINTYPE::iterator::operator->()->second); }
    BinF_t operator*() { return ANABINTYPE::iterator::operator*().second; };
  };
  ANABINITER begin() { return cont_.begin(); };
  ANABINITER end() { return cont_.end(); };
  class ANABINCITER : public ANABINTYPE::const_iterator
  {
   public:
    ANABINCITER() : ANABINTYPE::const_iterator() {};
    ANABINCITER(ANABINTYPE::const_iterator s) : ANABINTYPE::const_iterator(s) {};
    const BinF_t* operator->() const { return (BinF_t* const)&(ANABINTYPE::const_iterator::operator->()->second); }
    const BinF_t operator*() const { return ANABINTYPE::const_iterator::operator*().second; };
  };
  ANABINCITER begin() const { return cont_.begin(); };
  ANABINCITER end() const { return cont_.end(); };
  bool operator<(const AnaBin_t& a) const { return cont_ < a.cont_; }
  bool operator==(const AnaBin_t& a) const { return cont_ == a.cont_; }
   
  AnaBin_t() : cont_() {};
  AnaBin_t(const AnaBin_t& a) { setbin(a); };
  AnaBin_t(const ANABINTYPE& a) { setbin(a); };
  AnaBin_t(const BinF_t& a) { setbin(a); };
  AnaBin_t(const std::vector<AnaBin_t>& a) { setbin(a); };
  AnaBin_t(const std::string& n, const float& a, const float& b, const float& c=-99., const float& d=-99.) { setbin(n, a, b, c, d); };
  const BinF_t getbin (const std::string& name) const {
    const auto& b = cont_.find(name);
    if (b!=cont_.end()) { return b->second; }
    return BinF_t();
  };
  const BinF_t getbin (const int& i) const {
    return std::next(cont_.begin(), i)->second;
  };
  void setbin (const BinF_t& a, const bool& replace=true) {
    if (a.valid() && (replace || !hasbin(a.name()))) {
      cont_[a.name()] = a;
    }
  };
  void setbin (const ANABINTYPE& a, const bool& replace=true) {
    for (const auto& c : a) {
      setbin(c.second, replace);
    }
  };
  void setbin (const AnaBin_t& a, const bool& replace=true) {
    setbin(a.cont_, replace);
  };
  void setbin (const std::vector<AnaBin_t>& v, const bool& replace=true) {
    for (const auto& a : v) {
      setbin(a, replace);
    }
  };
  void setbin (const std::vector<BinF_t>& v, const bool& replace=true) {
    for (const auto& a : v) {
      setbin(a, replace);
    }
  };
  void setbin (const std::string& n, const float& a, const float& b, const float& c=-99., const float& d=-99.) {
    setbin(BinF_t(n,a,b,c,d));
  };
  void removebin (const std::string& name) {
    const auto& b = cont_.find(name);
    if (b!=cont_.end()) { cont_.erase(b); }
  };
  bool hasbin (const std::string& name) const {
    return (cont_.find(name)!=cont_.end());
  };
  bool hasbin (const BinF_t& bin) const {
    const auto& b = cont_.find(bin.name());
    return (b!=cont_.end() ? b->second==bin : false);
  };
  void print() const {
    std::string l = "";
    for (const auto& e : cont_) {
      const auto& s = e.second;
      l += (s.name()+" = [ "+std::to_string(s.mean())+" , "+std::to_string(s.width())+" , "+std::to_string(s.low())+" , "+std::to_string(s.high())+" ] , ");
    }
    if (l.rfind("] , ")!=std::string::npos) { l.erase(l.rfind(" , "), 3); }
    std::cout << l << std::endl;
  }
  size_t size() const { return cont_.size(); }
 private:
  ANABINTYPE cont_;
};


typedef std::pair<std::string, std::vector<float>  > FloatVecPair_t;
typedef std::vector< FloatVecPair_t                > FloatVecPairVec_t;
typedef std::map< AnaBin_t , FloatVecPairVec_t       > BinInfoMap_t;
typedef std::map< AnaBin_t , BinInfoMap_t            > BinInfoDiMap_t;
typedef std::vector< BinInfoDiMap_t                > BinInfoDiMapVec_t;
typedef std::map< std::string , BinInfoDiMapVec_t  > BinInfoTriMap_t;


#endif // #ifndef BIN_H_
