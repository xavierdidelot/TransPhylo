#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the probability of a phylogenetic tree given a transmission tree
//' @param ctree Combined phylogenetic/transmission tree
//' @param neg Within-host coalescent rate
//' @param w Vector of hosts for which to calculate the probability, or nothing for all
//' @return Probability of phylogeny given transmission tree
//' @export
// [[Rcpp::export]]
double probPTreeGivenTTree(NumericMatrix ctree, double neg, IntegerVector w = IntegerVector(0))  {
  double prob=0;
  double n=max(ctree(_,3));
  if (w.length()==0) w=seq(1,n);
  IntegerVector table=IntegerVector(n+1);
  IntegerVector parents=IntegerVector(ctree.nrow()+1);
  for (int j=0;j<ctree.nrow();j++) {
    table(ctree(j,3))++;
    parents(ctree(j,1))=j+1;
    parents(ctree(j,2))=j+1;
  }
  Environment pkg = Environment::namespace_env("TransPhylo");
  Function probSubtree = pkg["probSubtree"];
  for (int i=0;i<w.length();i++) {
    if (table[w(i)]<2) continue;
    std::vector<int> toinc;
    for (int j=0;j<ctree.nrow();j++) if (ctree(j,3)==w(i)) toinc.push_back(j+1);
    toinc.push_back(parents(toinc.back()));
    IntegerVector revtoinc=IntegerVector(ctree.nrow()+1);
    for (int j=0;j<toinc.size();j++) revtoinc[toinc[j]]=j+1;
    NumericMatrix subtree=NumericMatrix(toinc.size(),2);
    for (int j=0;j<toinc.size();j++) {
      subtree(j,0)=ctree(toinc[j]-1,0);
      subtree(j,1)=revtoinc(parents(toinc[j]));
    }
    RObject res = probSubtree(subtree,neg);
    prob += as<double>(res);
  } 
  return(prob);
}

// [[Rcpp::export]]
double coalescent(NumericVector leaves, NumericVector nodes, double alpha) {
  int n=leaves.length();
  double p = -log(alpha) * (n - 1);
  int i1=1,i2=0,k=1;
  double prev=leaves[0];
  for (int i=1;i<(n+n-1);i++) {
    if (i1<n&&leaves[i1]>nodes[i2]) {
      p = p - (k * (k - 1.0) / (2.0 * alpha) * (prev-leaves[i1]));
      prev=leaves[i1++];
      k++;
    } else {
      p = p - (k * (k - 1.0) / (2.0 * alpha) * (prev-nodes[i2]));
      prev=nodes[i2++];
      k--;
    }
  }
  return(p);
}

