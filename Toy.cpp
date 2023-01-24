#include "headToy.h"

int32 seed = (int32)time(0);

CRandomMersenne rg(seed);
	Pars pa;//parameter in struct
	initials in; //initial conditions
	totals to;
	Times ti;
	vector<Patch> Site; // information on each population
	vector<PatchGroup> Group; // information on each population
	ostringstream os;ofstream run;
	std::clock_t star;
	//PatchArray pat;						
	//VilArray vil;				
	//double inga;
	//double xstart,ystart,MaxDis;
	int main(void)
		{
			cin>>pa.set;
			cin>>in.inputfile;
			cin>>ti.N;//number of replicates
			cin>>ti.interval;
			cin>>ti.maxT;
			cin>>ti.rec;
			cin>>in.W;
			cin>>in.H;
			cin>>in.heg_time;// when gene drive is released
			cin>>pa.m;// migration rate
			cin>>pa.delta;// gene drive migration rate
			cin>>pa.b;// gene drive conversion of wildtype rate
			cin>>pa.xw;// wildtype extinction rate
			cin>>pa.xh;// gene drive extinction rate
			cin>>pa.LD;
			cin>>pa.ngroups;
			cin>>pa.U;
			star = std::clock();
			RunNReps(ti.N);
		return 0;};


	void RunNReps(int N)
	      {
		for(int i=0;i<N;i++)
		{
			initiate();
			RunMaxT(i);
			//EstRStar();			
		};
		return;};

	void RunMaxT(int i)
		{
		int TT=0;
		double TTT=0;
		double TTTT=0;
		int TTTTT=0;
		int swit=0;
		while (TT<ti.maxT/ti.interval)
			{
int th=0;
			TT++;
			TTT+=ti.interval;
			TTTT+=ti.interval;
			cout<<TTT-in.heg_time<<"      "<<to.W<<"       "<<to.H<<"       "<<to.E<<"       "<<to.NumPat<<"       "<<to.GsiI<<"     "<<to.GseS<<"     "<<to.GieI<<endl;
			if(TTT>=in.heg_time)
			{
				if(swit==0)
					{
					for(int index=0;index<in.H;index++){if(to.W>0) PutH();};
					swit=1;
					};
			//if(swit==1)cout<<TT*ti.interval-in.heg_time<<"       "<<to.H<<endl;
				if(TTTT>=ti.rec)
					{
					record(TTTTT);
					TTTTT++;
					TTTT-=ti.rec;
					};
			};
		if(TTT<in.heg_time ||(TTT>in.heg_time &&to.H>0))	RunOnceInt(ti.interval);
			};
        return;};



	void RunOnceInt(double interval)
	  {
	    double T;
	    T=0;
	    double dt;
	    while(T<interval)
		    {
		if(to.H+to.GseS>0)
			{
			dt=OneStep();
			T+=dt;
			if(to.GsiI<0.0000001)to.GsiI=0;
			if(to.GseS<0.0000001)to.GseS=0;
			if(to.GieI<0.0000001)to.GieI=0;
		//	cout<<" roi     "<<to.W<<"       "<<to.H<<"       "<<to.E<<"       "<<to.NumPat<<"       "<<to.GsiI<<"     "<<to.GseS<<"     "<<to.GieI<<endl;
			}
		else T=interval;
		    };
	return;}

	void initiate(){
		to.W=0;to.H=0;to.E=0;to.GseS=0;to.GsiI=0;to.NumPat=0;to.GieI=0;
		/*----------------------------------------------------------------------------------*/
		//cout<<"in1"<<endl;
		/*------------------------------input the settlement data---------------------------*/
		ostringstream ddd;
		ddd.str(in.inputfile);
		ifstream indata(ddd.str().c_str()); 
		string line ="";
		string cell,name,type;
		int site; double weight;
		Patch pp;
		PatchGroup gg;
		while( getline(indata, line ))
		{       
			stringstream lineStream(line);
			getline(lineStream,cell,','); pp.x=strtod(cell.c_str(),NULL);
			getline(lineStream,cell,','); pp.y=strtod(cell.c_str(),NULL);
			pp.GseS=0;
			pp.GsiI=0;
			pp.GieI=0;
			pp.type='E';
			pp.connecIND.clear();
			pp.connecW.clear();
			Site.push_back(pp);
			to.NumPat++;
			to.E++;
			//SetsPerCell[pp.group.][pp.group.].push_back(Site.size()-1);
		};
		/*----------------------------------------------------------------------------------*/
		
		for(int ii=0;ii<pa.ngroups;ii++)
		{
			gg.NumPat=0;
			gg.GsiI=0;
			gg.GseS=0;
			gg.GieI=0;
			gg.W=0;
			gg.H=0;
			gg.E=0;
			gg.members.clear();
			Group.push_back(gg);
		};
		int gp;
		for(int ii=0;ii<Site.size();ii++)
			{
				gp=rg.IRandomX(0,pa.ngroups-1);
				Site[ii].group=gp;
				Group[gp].NumPat++;
				Group[gp].E++;
				Site[ii].groupindex=Group[gp].NumPat-1;
				Group[gp].members.push_back(ii);
			};
		UpdateConnec();
		
		for(int index=0;index<in.W;index++) PutW();
//		for(int index=0;index<in.H;index++) PutH();
	return;};




	void UpdateConnec()
	{
	double dd,ww,cor;
	for(int index=0;index<Site.size();index++)
	{
		Site[index].connecIND.clear();
		Site[index].connecW.clear();
		Site[index].TotW=0;
		for(int ii=0;ii<Site.size();ii++)
			{
			if(ii!=index)
				{
				//dd=dist(Site[index].x,Site[index].y,Site[ii].x,Site[ii].y);
				dd=distB(pa.U,Site[index].x,Site[index].y,Site[ii].x,Site[ii].y);
				if(dd<pa.LD)
				{
					Site[index].connecIND.push_back(ii); 
					ww=(1-dd/pa.LD)*(3/(PI*pa.LD*pa.LD))*(pa.U*pa.U/to.NumPat);
					Site[index].connecW.push_back(ww); 
					Site[index].TotW+=ww;
//					cout<<index<<"   "<<ii<<endl;
				};
				};
		
			};
	};

	return;};

	void PutW()
	{
		int ind,indB,indC,gp;
		double ww;
		int test=0;
		while(test==0)
			{
			gp=rg.IRandomX(0,pa.ngroups-1);
			if(double(1.0*Group[gp].E/(1.0*to.E))>rg.Random())test=1;
			};
		test=0;
		while(test==0)
			{
			ind=rg.IRandomX(0,Group[gp].NumPat-1);
			if(Site[Group[gp].members[ind]].type=='E')test=1;
			};
	
		indB=Group[gp].members[ind];
		Site[indB].type='W';
		int howmany=Site[indB].connecIND.size();
//		cout<<"pw "<<indB<<"   "<<gp<<"   "<<howmany<<endl;
		for(int ii=0;ii<howmany;ii++)
		{	
		indC=Site[indB].connecIND[ii];	
		ww=Site[indB].connecW[ii];	

		if(Site[indC].type=='E')
			{
			Site[indC].GseS+=ww;
			Group[Site[indC].group].GseS+=ww;
			to.GseS+=ww;
		//cout<<"pw "<<indB<<"   "<<indC<<"   "<<ww<<endl;
			};
		if(Site[indC].type=='H')
			{
			Site[indB].GsiI+=ww;
			Group[gp].GsiI+=ww;
			to.GsiI+=ww;
			};
		};
		to.W++;
		Group[gp].W++;
		to.E--;
		Group[gp].E--;
		Group[gp].GseS-=Site[indB].GseS;
		to.GseS-=Site[indB].GseS;
		Group[gp].GieI-=Site[indB].GieI;
		to.GieI-=Site[indB].GieI;
		Site[indB].GseS=0;
		Site[indB].GieI=0;
	return;};

	void PutH()
	{
		int test=0;
		int ind,indB,indC,gp;
		double ww;
		//cout<<"ph1"<<endl;
		while(test==0)
			{
			gp=rg.IRandomX(0,pa.ngroups-1);
			if(double((1.0*Group[gp].E+Group[gp].W)/(1.0*to.NumPat))>rg.Random())test=1;
			};
		test=0;
		while(test==0)
			{
			ind=rg.IRandomX(0,Group[gp].NumPat-1);
			if(Site[Group[gp].members[ind]].type=='E' || Site[Group[gp].members[ind]].type=='W')test=1;
			};

		indB=Group[gp].members[ind];
		int howmany=Site[indB].connecIND.size();
		for(int ii=0;ii<howmany;ii++)
		{	
		indC=Site[indB].connecIND[ii];	
		ww=Site[indB].connecW[ii];	

		if(Site[indC].type=='E' && Site[indB].type=='W')
			{
			Site[indC].GseS-=ww;
			Group[Site[indC].group].GseS-=ww;
			to.GseS-=ww;
			};
		if(Site[indC].type=='E')
			{
			Site[indC].GieI+=ww;
			Group[Site[indC].group].GieI+=ww;
			to.GieI+=ww;
			};
		if(Site[indC].type=='W')
			{
			Site[indC].GsiI+=ww;
			Group[Site[indC].group].GsiI+=ww;
			to.GsiI+=ww;
			};
		};

		if(Site[indB].type=='W')
			{
			Group[gp].W--;
			to.W--;
			Group[gp].GsiI-=Site[indB].GsiI;
			to.GsiI-=Site[indB].GsiI;
			Site[indB].GsiI=0;
			}
		else
			{
			to.E--;
			Group[gp].E--;
			Group[gp].GseS-=Site[indB].GseS;
			to.GseS-=Site[indB].GseS;
			Group[gp].GieI-=Site[indB].GieI;
			to.GieI-=Site[indB].GieI;
			Site[indB].GseS=0;
			Site[indB].GieI=0;
			};
		Site[indB].type='H';
		Group[gp].H++;
		to.H++;
	      return;};



	double OneStep()
	{
		//	cout<<" os1   "<<to.W<<"       "<<to.H<<"       "<<to.E<<"       "<<to.NumPat<<"       "<<to.GsiI<<"     "<<to.GseS<<"     "<<to.GieI<<endl;
		double rate[4];
		rate[0] = pa.b*pa.m*to.GsiI;//W to H
		rate[1] = to.W*pa.xw;//W to E
		rate[2] = to.H*pa.xh;//H to E
		rate[3] = pa.m*to.GseS;//E to W
		rate[4] = pa.m*pa.delta*to.GieI;//E to H
		double sum=rate[0]+rate[1]+rate[2]+rate[3]+rate[4];
		//double dt=ti.interval;
		double dt=ti.interval;
		if(sum>0)
		{
		dt=random_exp(sum);
	//	cout<<rate[0]<<"   "<<rate[1]<<"   "<<rate[2]<<"   "<<rate[3]<<"   "<<rate[4]<<endl;
			int i= pick(rate,4,sum);
	//			cout<<" os1     "<<i<<"  "<<to.W<<"       "<<to.H<<"       "<<to.E<<"       "<<to.NumPat<<"       "<<to.GsiI<<"     "<<to.GseS<<"     "<<to.GieI<<endl;
			if(i==0) WtoH();
			else if(i==1) WtoE();
			else if(i==2) HtoE();
			else if(i==3)EtoW();
			else EtoH();
			if(to.W<0 || to.H<0 ||to.E<0)
				{cout<<" os2     "<<i<<"  "<<to.W<<"       "<<to.H<<"       "<<to.E<<"       "<<to.NumPat<<"   gsi    "<<to.GsiI<<"  gse   "<<to.GseS<<"  gie   "<<to.GieI<<endl;
				exit(1);};
		};
	return dt;};



	void WtoH()
	  {
		int test=0;
		int ind,indB,indC,gp;
		double ww;
		while(test==0)
			{
			gp=rg.IRandomX(0,pa.ngroups-1);
			if(Group[gp].GsiI/to.GsiI>rg.Random())test=1;
			};

		test=0;
		while(test==0)
			{
			ind=rg.IRandomX(0,Group[gp].NumPat-1);
			if(Site[Group[gp].members[ind]].GsiI/Group[gp].GsiI>rg.Random())test=1;
			};

		indB=Group[gp].members[ind];
		int howmany=Site[indB].connecIND.size();
		for(int ii=0;ii<howmany;ii++)
		{	
		indC=Site[indB].connecIND[ii];	
		ww=Site[indB].connecW[ii];	

		if(Site[indC].type=='E')
			{
			Site[indC].GseS-=ww;
			Group[Site[indC].group].GseS-=ww;
			to.GseS-=ww;
			Site[indC].GieI+=ww;
			Group[Site[indC].group].GieI+=ww;
			to.GieI+=ww;
			};
		if(Site[indC].type=='W')
			{
			Site[indC].GsiI+=ww;
			Group[Site[indC].group].GsiI+=ww;
			to.GsiI+=ww;
			};
		};

		Group[gp].W--;
		to.W--;
		Group[gp].GsiI-=Site[indB].GsiI;
		to.GsiI-=Site[indB].GsiI;
		Site[indB].GsiI=0;
		Site[indB].type='H';
		Group[gp].H++;
		to.H++;
	      return;};

	void HtoE()
	  {
		int test=0;
		int ind,indB,indC,gp;
		double ww;
		while(test==0)
			{
			gp=rg.IRandomX(0,pa.ngroups-1);
			if(Group[gp].H/(1.0*to.H)>rg.Random())test=1;
			};

		test=0;
		while(test==0)
			{
			ind=rg.IRandomX(0,Group[gp].NumPat-1);
			if(Site[Group[gp].members[ind]].type=='H')test=1;
			};

		indB=Group[gp].members[ind];
		int howmany=Site[indB].connecIND.size();
		for(int ii=0;ii<howmany;ii++)
		{	
		indC=Site[indB].connecIND[ii];	
		ww=Site[indB].connecW[ii];	
		if(Site[indC].type=='W')
			{
			Site[indB].GseS+=ww;
			Group[gp].GseS+=ww;
			to.GseS+=ww;
			Site[indC].GsiI-=ww;
			Group[Site[indC].group].GsiI-=ww;
			to.GsiI-=ww;
			};
		if(Site[indC].type=='E')
			{
			Site[indC].GieI-=ww;
			Group[Site[indC].group].GieI-=ww;
			to.GieI-=ww;
			};
		if(Site[indC].type=='H')
			{
			Site[indB].GieI+=ww;
			Group[gp].GieI+=ww;
			to.GieI+=ww;
			};
		};
		Site[indB].type='E';
		Group[gp].H--;
		to.H--;
		Group[gp].E++;
		to.E++;
	      return;};

	void WtoE()
	  {

		int test=0;
		int ind,indB,indC,gp;
		double ww;
		while(test==0)
			{
			gp=rg.IRandomX(0,pa.ngroups-1);
			if(Group[gp].W/(1.0*to.W)>rg.Random())test=1;
			};
		test=0;
		while(test==0)
			{
			ind=rg.IRandomX(0,Group[gp].NumPat-1);
			if(Site[Group[gp].members[ind]].type=='W')test=1;
			};

		indB=Group[gp].members[ind];
		int howmany=Site[indB].connecIND.size();
		for(int ii=0;ii<howmany;ii++)
		{	
		indC=Site[indB].connecIND[ii];	
		ww=Site[indB].connecW[ii];	
		if(Site[indC].type=='E')
			{
			Site[indC].GseS-=ww;
			Group[Site[indC].group].GseS-=ww;
			to.GseS-=ww;
			};
		if(Site[indC].type=='W')
			{
			Site[indB].GseS+=ww;
			Group[gp].GseS+=ww;
			to.GseS+=ww;
			};
		if(Site[indC].type=='H')
			{
			Site[indB].GieI+=ww;
			Group[gp].GieI+=ww;
			to.GieI+=ww;
			};
		};

		Group[gp].W--;
		to.W--;
		Group[gp].GsiI-=Site[indB].GsiI;
		to.GsiI-=Site[indB].GsiI;
		Site[indB].GsiI=0;
		Site[indB].type='E';
		Group[gp].E++;
		to.E++;
	      return;};

	void EtoW()
	  {
		int test=0;
		int ind,indB,indC,gp;
		double ww;
		while(test==0)
			{
			gp=rg.IRandomX(0,pa.ngroups-1);
			if(Group[gp].GseS/to.GseS>rg.Random())test=1;
			};
		test=0;
		while(test==0)
			{
			ind=rg.IRandomX(0,Group[gp].NumPat-1);
			if(Site[Group[gp].members[ind]].GseS/Group[gp].GseS>rg.Random())test=1;
			};

		indB=Group[gp].members[ind];
		int howmany=Site[indB].connecIND.size();
		for(int ii=0;ii<howmany;ii++)
		{	
		indC=Site[indB].connecIND[ii];	
		ww=Site[indB].connecW[ii];	
		if(Site[indC].type=='E')
			{
			Site[indC].GseS+=ww;
			Group[Site[indC].group].GseS+=ww;
			to.GseS+=ww;
			};
		if(Site[indC].type=='H')
			{
			Site[indB].GsiI+=ww;
			Group[gp].GsiI+=ww;
			to.GsiI+=ww;
			};
		};
		to.W++;
		Group[gp].W++;
		to.E--;
		Group[gp].E--;
		Group[gp].GseS-=Site[indB].GseS;
		to.GseS-=Site[indB].GseS;
		Site[indB].GseS=0;
		Group[gp].GieI-=Site[indB].GieI;
		to.GieI-=Site[indB].GieI;
		Site[indB].GieI=0;
		Site[indB].type='W';

	      return;};


	void EtoH()
	  {
		int test=0;
		int ind,indB,indC,gp;
		double ww;
		while(test==0)
			{
			gp=rg.IRandomX(0,pa.ngroups-1);
			if(Group[gp].GieI/to.GieI>rg.Random())test=1;
			};
		test=0;
		while(test==0)
			{
			ind=rg.IRandomX(0,Group[gp].NumPat-1);
			if(Site[Group[gp].members[ind]].GieI/Group[gp].GieI>rg.Random())test=1;
			};

		indB=Group[gp].members[ind];
		int howmany=Site[indB].connecIND.size();
		for(int ii=0;ii<howmany;ii++)
		{	
		indC=Site[indB].connecIND[ii];	
		ww=Site[indB].connecW[ii];	

		if(Site[indC].type=='E')
			{
			Site[indC].GieI+=ww;
			Group[Site[indC].group].GieI+=ww;
			to.GieI+=ww;
			};
		if(Site[indC].type=='W')
			{
			Site[indC].GsiI+=ww;
			Group[Site[indC].group].GsiI+=ww;
			to.GsiI+=ww;
			};
		};

		to.E--;
		Group[gp].E--;
		Group[gp].GseS-=Site[indB].GseS;
		to.GseS-=Site[indB].GseS;
		Site[indB].GseS=0;
		Group[gp].GieI-=Site[indB].GieI;
		to.GieI-=Site[indB].GieI;
		Site[indB].GieI=0;
		Site[indB].type='H';
		Group[gp].H++;
		to.H++;

	      return;};


	int pick (int* list,int len,int sum)
	{
		int test=0;
		int i;
		while(test==0)
			{i=rg.IRandomX(0,len);
			if(rg.Random()<list[i]/sum)test=1;};
	return i;};

	int pick (double* list,int len,double sum)
	{
		int test=0;
		int i;
		while(test==0)
			{i=rg.IRandomX(0,len);
			if(rg.Random()<list[i]/sum)test=1;};
	return i;};


double dmodulo(double i, double n)
	{ return fmod((fmod(i,n)+n),n); };

int random_poisson(double landa)
	{
	double p=exp(-landa);
	double g=p;
	double u=rg.Random();
	int k=0;
	while (u>g)
	    {
	    p*=(landa/(double)(++k));
		g+=p;
	    };
	return k;
	};



double random_exp(double lamda)
{double dd= rg.Random()*0.9999;
return (-log(1.0-dd)/(double)lamda);
};

double dist (double x1, double y1, double x2, double y2)
	{
		return double(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
	};


double distB (double U,double x1, double y1, double x2, double y2)
	{
		double xdist=0;
		double ydist=0;

		if(abso(x1-x2)>U-abso(x1-x2)) xdist=U-abso(x1-x2);
		if(abso(x1-x2)<U-abso(x1-x2)) xdist=abso(x1-x2);
		if(abso(y1-y2)>U-abso(y1-y2)) ydist=U-abso(y1-y2);
		if(abso(y1-y2)<U-abso(y1-y2)) ydist=abso(y1-y2);

		return double(sqrt(xdist*xdist+ydist*ydist));
	};

double abso(double XX)
{double YY=XX;
	if(XX<0)YY*=-1;
	return YY;};


void record(int T)
{
	ostringstream os;
	ofstream logfile;
	os<<"PRes"<<pa.set<<"T"<<T<<".txt";
	logfile.open(os.str().c_str());

	for(int i=0;i<Site.size();i++)
					{
			//logfile<<Site[i].x<<"     "<<Site[i].y<<"     "<<Site[i].type<<endl;
			logfile<<Site[i].type<<endl;
					};
	os.str("");
	logfile.close();
	return;};



	void CRandomMersenne::Init0(uint32 seed) {
   // Detect computer architecture
   union {float f; uint32 i[2];} convert;
   convert.f = 1.0;
   if (convert.i[1] == 0x3FF00000) Architecture = LITTLE_ENDIAN1;
   else if (convert.i[0] == 0x3FF00000) Architecture = BIG_ENDIAN1;
   else Architecture = NONIEEE;

   // Seed generator
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(uint32 seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int i = 0; i < 37; i++) BRandom();
}


void CRandomMersenne::RandomInitByArray(uint32 seeds[], int length) {
   // Seed by more than 32 bits
   int i, j, k;

   // Initialize
   Init0(19650218);

   if (length <= 0) return;

   // Randomize mt[] using whole seeds[] array
   i = 1;  j = 0;
   k = (MERS_N > length ? MERS_N : length);
   for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + seeds[j] + j;
      i++; j++;
      if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
      if (j >= length) j=0;}
   for (k = MERS_N-1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
      if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
   mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

   // Randomize some more
   mti = 0;
   for (int i = 0; i <= MERS_N; i++) BRandom();
}


uint32 CRandomMersenne::BRandom() {
   // Generate 32 random bits
   uint32 y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32 LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32 UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32 mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }

   y = mt[mti++];

#if 1
   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;
#endif

   return y;
}


double CRandomMersenne::Random() {
   // Output random float number in the interval 0 <= x < 1
   union {float f; uint32 i[2];} convert;
   double ra=1.1;
	while(ra>=1)
	{
   uint32 r = BRandom();               // Get 32 random bits
   // The fastest way to convert random bits to floating point is as follows:
   // Set the binary exponent of a floating point number to 1+bias and set
   // the mantissa to random bits. This will give a random number in the
   // interval [1,2). Then subtract 1.0 to get a random number in the interval
   // [0,1). This procedure requires that we know how floating point numbers
   // are stored. The storing method is tested in function RandomInit and saved
   // in the variable Architecture.

   // This shortcut allows the compiler to optimize away the following switch
   // statement for the most common architectures:
#if defined(_M_IX86) || defined(_M_X64) || defined(__LITTLE_ENDIAN__)
   Architecture = LITTLE_ENDIAN1;
#elif defined(__BIG_ENDIAN__)
   Architecture = BIG_ENDIAN1;
#endif

   switch (Architecture) {
   case LITTLE_ENDIAN1:
      convert.i[0] =  r << 20;
      convert.i[1] = (r >> 12) | 0x3FF00000;
      return convert.f - 1.0;
   case BIG_ENDIAN1:
      convert.i[1] =  r << 20;
      convert.i[0] = (r >> 12) | 0x3FF00000;
      return convert.f - 1.0;
   case NONIEEE: default: ;
   }
   // This somewhat slower method works for all architectures, including
   // non-IEEE floating point representation:
   //return 0.0000005 +0.999999*(double)r * (1./((double)(uint32)(-1L)+1.));
   ra= (double)r * (1./((double)(uint32)(-1L)+1.));
	};
   return ra;
}


int CRandomMersenne::IRandom(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Multiply interval with random and truncate
   int r = int((max - min + 1) * Random()) + min;
   if (r > max) r = max;
   return r;
}


int CRandomMersenne::IRandomX(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Each output value has exactly the same probability.
   // This is obtained by rejecting certain bit values so that the number
   // of possible bit values is divisible by the interval length
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
#ifdef  INT64_DEFINED
   // 64 bit integers available. Use multiply and shift method
   uint32 interval;                    // Length of interval
   uint64 longran;                     // Random bits * interval
   uint32 iran;                        // Longran / 2^32
   uint32 remainder;                   // Longran % 2^32

   interval = uint32(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder = 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then
      RLimit = uint32(((uint64)1 << 32) / interval) * interval - 1;
      LastInterval = interval;
   }
   do { // Rejection loop
      longran  = (uint64)BRandom() * interval;
      iran = (uint32)(longran >> 32);
      remainder = (uint32)longran;
   } while (remainder > RLimit);
   // Convert back to signed and return result
   return (int32)iran + min;

#else
   // 64 bit integers not available. Use modulo method
   uint32 interval;                    // Length of interval
   uint32 bran;                        // Random bits
   uint32 iran;                        // bran / interval
   uint32 remainder;                   // bran % interval

   interval = uint32(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when iran = 2^32 / interval
      // We can't make 2^32 so we use 2^32-1 and correct afterwards
      RLimit = (uint32)0xFFFFFFFF / interval;
      if ((uint32)0xFFFFFFFF % interval == interval - 1) RLimit++;
   }
   do { // Rejection loop
      bran = BRandom();
      iran = bran / interval;
      remainder = bran % interval;
   } while (iran >= RLimit);
   // Convert back to signed and return result
   return (int32)remainder + min;

#endif
}
