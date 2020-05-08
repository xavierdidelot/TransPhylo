#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double probSubtree(NumericMatrix tab, double rate) {
//tab[,0]=times at bottom;tab[,1]=father;rate=coalescence rate 
//Return the log-prior probability of a subtree 
//This is an extension to Eq1 of Drummond et al(2002) Genetics 161:1307-1320 that accounts for condition TMRCA<INCUBATION_PERIOD 
 int i,k,ii,c,l,li,anc;
 double p=0;
 double ma=max(tab(_,0));
 NumericVector tab0=ma-tab(_,0);//convert times to ages
 NumericVector tab1=tab(_,1)-1;//index start at zero
 LogicalVector isiso(tab.nrow(),true);
 for (i=0;i<tab.nrow()-1;i++) isiso[tab1[i]]=false;
 int niso=sum(isiso);
 IntegerVector ex(tab.nrow(),0);
 IntegerVector ind(tab.nrow(),0);
 c=0;
 for (i=1;i<tab.nrow();i++) {
   if (tab0[i]<tab0[0]) ind[c]=i; else ind [1+c]=i;
   c++;
 }
 IntegerVector iso(niso,0);
 c=0;
 for (i=1;i<tab.nrow();i++) {
   if (isiso[i]) {if (tab0[i]<tab0[0]) iso[c]=i; else iso[1+c]=i;
   c++;}
 }
 int cur=iso[0];//Start with youngest leaf
 LogicalVector isanc(tab.nrow(),false);//Ancestors of the current node
 while (tab1[cur]>=0) {ex[cur]=1;cur=tab1[cur];}//Activate path to root
 ex[ex.length()-1]=1;//Activate root
 bool start,found;
 double bra1,bra2,curage;
 for (li=1;li<iso.length();li++) {//For all leaves in increasing order of age
   l=iso[li];
   isanc.fill(false);
   anc=l;
   while (tab1[anc]>=0) {
     isanc[anc]=true;
     anc=tab1[anc];
   }
   bra1=0;
   bra2=0;
   start=false;
   found=false;
   curage=0;
   k=0;
   for (ii=0;ii<ind.length();ii++) {
     i=ind[ii];
     if (i==l) {start=true;curage=tab0[l];continue;}
     if (ex[i]==0) continue;//Ignore non-existent nodes
     if (start) {
       if (!found) {
         bra1+=k*(tab0[i]-curage);
         if (isanc[i]) found=true;
       }
       bra2+=k*(tab0[i]-curage);
     }
     curage=tab0[i];
     if (isiso[i]) k++; else k=k-ex[i]+1;
   }
   p=p-log(rate)-bra1/rate-log(1-exp(-bra2/rate));
   //if (k != 1) Rprintf("errorHere %d %f %f\n",k,bra1,bra2);
   //if (bra1 >= bra2)  Rprintf("errorThere %f %f\n",bra1,bra2);
   
   //Make all ancestors of current node active
   cur=l;
   while (1) {
    if (ex[cur]==1) {ex[cur]=2;break;}
     ex[cur]=1;
    cur=tab1[cur];
   }
  }
   return(p);
}

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
  //Environment pkg = Environment::namespace_env("TransPhylo");
  //Function probSubtree = pkg["probSubtree"];
  for (int i=0;i<w.length();i++) {
    if (table[w(i)]<2) continue;
    std::vector<int> toinc;
    for (int j=0;j<ctree.nrow();j++) if (ctree(j,3)==w(i)) toinc.push_back(j+1);
    toinc.push_back(parents(toinc.back()));
    IntegerVector revtoinc=IntegerVector(ctree.nrow()+1);
    for (unsigned int j=0;j<toinc.size();j++) revtoinc[toinc[j]]=j+1;
    NumericMatrix subtree=NumericMatrix(toinc.size(),2);
    for (unsigned int j=0;j<toinc.size();j++) {
      subtree(j,0)=ctree(toinc[j]-1,0);
      subtree(j,1)=revtoinc(parents(toinc[j]));
    }
    //RObject res = probSubtree(subtree,neg);
    //prob += as<double>(res);
    prob += probSubtree(subtree,neg);
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

