#ifndef Utilities_bin_h
#define Utilities_bin_h

#include <iostream>
#include <tuple>
#include <string>
#include <set>

// a simple template class to store a bin and overload the equality operator

// define a few common uses of the template class
template <typename T> class bin : public std::tuple<std::string,T,T>
{
 public:
  bin() : std::tuple<std::string,T,T>() {};
  bin(std::string n, T a, T b) : std::tuple<std::string,T,T>(n,a,b) {};
  std::string name() const { return std::get<0>(*this); }
  T low()  const { return std::get<1>(*this); }
  T high() const { return std::get<2>(*this); }
};
typedef bin<float>  binF;

// associate a set of bins to make an analysis bin
class anabin : public std::set<binF> {
   public:
      anabin() : std::set<binF>() {};
      anabin(const std::set<binF>& set) : std::set<binF>(set) {};
      binF getbin (const std::string& name) const {
	for (const auto& s : *this) { if (s.name()==name) { return s; } };
	return binF("", -1., -1.);
      };
      void setbin (const binF& bin) { this->insert(bin); };
      void setbin (const std::string& n, const float& a, const float& b) { this->insert(binF(n,a,b)); };
      void print() const {
	std::string l = "";
	for (const auto& s : *this) { l += (s.name()+" = [ "+std::to_string(s.low())+" , "+std::to_string(s.high())+" ] , "); }
	if (l.rfind("] , ")!=std::string::npos) { l.erase(l.rfind(" , "), 3); }
	std::cout << l << std::endl;
      }
};

#endif // #ifndef bin_h
