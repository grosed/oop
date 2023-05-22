
#include <iostream>
#include <tuple>
#include <functional>
#include <algorithm>
#include <ranges>
#include <set>
#include <numeric>

#include <concepts>
#include <vector>

#include <fstream>
#include <sstream>
#include <string>

template <typename T>
T lexical_cast(const std::string& str)
{
    T var;
    std::istringstream iss;
    iss.str(str);
    iss >> var;
    // deal with any error bits that may have been set on the stream
    return var;
}

template<typename data_type,
template <typename... table_type_args> typename table_type,
template <typename... row_type_args> typename row_type>
table_type<row_type<data_type> > csvtable(const std::string& filename)
{
  table_type<row_type<data_type> > table;
  std::ifstream infile(filename);
  while(infile)
    {
      std::string s;
      if(!getline(infile,s)) break;
      std::istringstream ss(s);
      row_type<data_type> row;
      while(ss)
      {
         std::string s;
         if(!getline(ss,s,',')) break;
         row.push_back(lexical_cast<data_type>(s)); 
      }
      table.push_back(row);
    }
  return table;
}

#include <list>
#include <vector>

#define readcsv csvtable<double,std::vector,std::vector>


struct sumstats
{
  std::vector<double> SX;
  std::vector<double> SXX;
  sumstats(const std::vector<double>&);
  std::vector<double> getSX() const;
  std::vector<double> getSXX() const;
};

sumstats::sumstats(const std::vector<double>& X)
{
  SX = std::vector<double>(X.size()+1,0.0);
  SXX = std::vector<double>(X.size()+1,0.0);
  std::transform(X.begin(),X.end(),SX.begin()+1,[](const auto& x){return x;});
  std::partial_sum(SX.begin(), SX.end(),SX.begin());
  std::transform(X.begin(),X.end(),SXX.begin()+1,[](const auto& x){return x*x;});
  std::partial_sum(SXX.begin(), SXX.end(),SXX.begin());
}

std::vector<double> sumstats::getSX() const
{
  return SX;
}

std::vector<double> sumstats::getSXX() const
{
  return SXX;
}

struct normal_mean
{
  sumstats S;  
  normal_mean(const std::vector<double>&);
  double operator()(const int&,const int&);
};

normal_mean::normal_mean(const std::vector<double>& X) : S(X) {}

double normal_mean::operator()(const int& i,const int& j)
{
  double val = S.SX[j+1] - S.SX[i];
  val *= val;
  val /= (j - i + 1);
  val = -val;
  val += S.SXX[j+1] - S.SXX[i];
  return val;
}

template <typename Rtype>
auto split(std::invocable<const Rtype&,const Rtype&> auto f,
	   std::invocable<const Rtype&,const Rtype&> auto g,
	   const std::set<Rtype>& R)
{
  auto ra = *std::min_element(R.begin(),R.end());
  auto rb = *std::max_element(R.begin(),R.end());
  auto cost = [&](const Rtype& i){return std::pair(i,f(ra,i) + g(i,rb));};
  return std::ranges::min(R | std::views::transform(cost),[](const auto& a,const auto& b){return a.second < b.second;});
}

int main()
{

  /*
  X<int,double> g;
  std::set<int> R{0,1,2}; 
  auto res = split(g,g,R);
  */
  
  // std::vector<double> X {1,2,3,4};
  /*
  sumstats S(X);
  std::cout << "****************" << std::endl;
  for(const auto& val : S.getSX())
    {
      std::cout << val << std::endl;
    }
  std::cout << "****************" << std::endl;
    for(const auto& val : S.getSXX())
    {
      std::cout << val << std::endl;
    }
  */


  std::string filename {"X.csv"};
  // std::vector<std::vector<double> > X {readcsv(filename)};
  std::vector<std::vector<double> > table {readcsv(filename)};
  std::vector<double> X(table.size());
  for(int i = 0; i < table.size(); i++)
    {
      X[i] = table[i][0];
    }
  
  normal_mean cost(X);
  std::cout << cost(1,3) << std::endl;
  std::cout << cost(0,0) << std::endl;
  std::cout << cost(X.size()-1,X.size()-1) << std::endl;

  std::set<int> R;
  for(int i = 0; i < X.size(); i++)
    {
      R.insert(i);
    }
  auto res = split(cost,cost,R);
  std::cout << res.first << " : " << res.second << std::endl;
  
  return 0;
}
