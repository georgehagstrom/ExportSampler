// Here we make a hierarchical model of the sediment trap data in Mouw et. al.
//




// Functions Block

functions {



}

// Data
//
data {

int<lower=0> NObs; // Number of observations
int<lower=0> NGrids; // Number of grid boxes
int<lower=0> BiomeVec[NGrids];
int<lower=0> RegionVec[NGrids];

vector<lower=0>[NObs] POCFlux; //Measured POC Flux
vector<lower=0>[NObs] Depth; // Depth in Meters
vector[NGrids] POC114Model; // Modelled POC Flux at 114 meters depth


vector<lower=0>[NGrids] POC114GM; // Modelled POC Flux at 114 meters depth
vector<lower=0>[NGrids] POC114Redfield; // Modelled POC Flux at 114 meters depth

int<lower=0> IndVec[NObs]; // Grid location of each observation


}


transformed data {
vector[NObs] POC114ModelFull;

int NBiomes;
int NRegions;

NBiomes = 3;
NRegions = 13;

for (n in 1:NObs)
{
POC114ModelFull[n] = (POC114Model[IndVec[n]]);

}

}


parameters {


real <lower=0> sigmaPOC;
vector<lower=.4,upper=2.0>[3] kVec;
real<lower=0> sigmaK;
//real<lower=0.4,upper = 2.0> k;
real<lower=0,upper=10> sigmaMeasured;
vector<lower=0>[NObs] POCDepth;
}


transformed parameters {


vector<lower=0>[NObs] sigmaVec;
vector<lower=0>[NObs] DepthCorrection;
vector[NObs] POCReal;
real z0;
z0 = 114.0;

for (n in 1:NObs){

DepthCorrection[n] = 1.0/( (Depth[n]/z0)^(-kVec[BiomeVec[IndVec[n]]]));
sigmaVec[n] = sigmaPOC;
//DepthCorrection[n] = 1.0/( (Depth[n]/z0)^(-k));

}

POCReal = (POCDepth).*DepthCorrection;

}



model {


//sigmaMeasured ~ exponential(1.0);
sigmaPOC ~ exponential(20);
sigmaK ~ exponential(0.3);
kVec ~ normal(0.85,sigmaK);
//k ~ normal(0.85,0.3);


for (n in 1:NObs){
POCDepth[n] ~ normal(POCFlux[n],0.01) T[0,];
POC114ModelFull[n] ~ normal(POCReal[n],sigmaVec[n]);
}


}








generated quantities {





vector[NObs] log_lik;
vector[NObs] y_hat;

for (n in 1:NObs){

log_lik[n] = lognormal_lpdf(POC114ModelFull[n] | POCReal[n]   , sigmaVec[n]   );
y_hat[n] = lognormal_rng(POCReal[n],sigmaVec[n]);


} 





}



