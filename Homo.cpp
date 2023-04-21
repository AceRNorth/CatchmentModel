#include "head.h"

	/*--------------------------------------global variables----------------------------------*/
	/*----------------------------------------------------------------------------------------*/
	/*-------------random number seed-----------------*/
	int32 seed = (int32)time(0);
	CRandomMersenne rg(seed);
	/*------------------------------------------------*/
	/*----------------global structs------------------*/
	Pars pa;//parameters 
	initials in; //initial conditions
	totals to; // totals
	Times ti; // timekeeping parameters
	vector<Patch> Site; // information on each population
	/*------------------------------------------------*/
	/*------------------Output files------------------*/
	ostringstream os1,os2,os3,os4;
	ofstream globalinfo,localinfo,duration,ParList;
	/*------------------------------------------------*/
	/*------------Array to store rain input data--------------*/
//	double rain[nx][ny][52];//nx is x grid cell,ny is y grid cell, 52 is weeks in year
	/*--------------------------------------------------------*/

	/*----------------------------------------------------------------------------------------*/


	//	clock_t now=clock();
	int main(void)
		{
	/*---------input parameters---------------------*/
		cin>>pa.set; //enumeration for output files
		cin>>pa.index;
//		cin>>in.inputfile;// name of file with settlement data
		cin>>pa.NumPat;// number of sites
//		cin>>in.polyconnect;// name of file listing connected settlements per settlement
//		cin>>in.polyedges;// name of file listing connection weights
		cin>>in.recSitesFreq;//if outputting local data, how many sites to collect for (1 is all sites, 10 is 1 in 10 etc)
		cin>>ti.interval;//time interval for outputting global data
		cin>>ti.start;
		cin>>ti.maxT;// maximum time to simulate until (days)
		cin>>ti.recfreq;// how often to collect local data (every rec days)
		cin>>ti.recstart;// how often to collect local data (every rec days)
		cin>>ti.recend;// how often to collect local data (every rec days)
		cin>>ti.NumRuns;// how many simulation runs
		cin>>in.driver_time;// time to start releasing drive alleles; if negative, the release dates are random, if positive, the releases are on the specific date
		cin>>in.driver_start;// time to start releasing drive alleles
		cin>>in.driver_end;// time to end releasing drive alleles
		cin>>in.NumDriver;// number of drive homozygous male mosquitoes per release
		cin>>in.NumDriverSites;// number of release sites per year
		cin>>pa.muJ;//juvenile density independent mortality per day
		cin>>pa.muA;//adult mortality per day
		cin>>pa.d;//adult dispersal rate
		cin>>pa.gamma;//rate of r2 allele formation from W/D meiosis in females
		cin>>pa.Mgamma;//rate of r2 allele formation from W/D meiosis in males
		cin>>pa.beta;//parameter that controls mating rate
		cin>>pa.theta;//egg laying rate of wildtype females (eggs per day)
		cin>>pa.xi;// somatic Cas9 expression fitness cost
		cin>>pa.e;// homing rate in females
		cin>>pa.LD;//maximum distance at which settlements are connected

		cin>>pa.psi;//Aestivation rate
		cin>>pa.muAES;//aestivation mortality
		cin>>pa.t_hide1;//start day of going into aestivation
		cin>>pa.t_hide2;//end day of going into aestivation
		cin>>pa.t_wake1;//start day of emerging from aestivation
		cin>>pa.t_wake2;//end day of emerging from aestivation
		cin>>pa.alpha0;//baseline contribution to carrying capacity
		cin>>pa.alpha1;//Maximum contribution from rain
		cin>>pa.alpha2;//maximum contribution from water bodies
		cin>>pa.phi;// Increase in carrying capacity/alpha_1 per mm rain per week (when rainfall low)
		cin>>pa.delta;//Increase in length of standing water from non-permanent waterways per km non-permanent waterways (within 5km)per mm rain per week (when rainfall low)
		cin>>pa.kappa;//Increase in carrying capacity/α2 per km standing water (within 5km; when water bodies rare)
		cin>>pa.al0var;// variance in baseline carrying capacity
		cin>>pa.omega;//human density at which further increases do not increase carrying capacity
//		cin>>in.rainfile;// rain data input file
//		cin>>in.mortfile;// mortality data input file
		cin>>pa.bias;// male proportion in progeney of males carrying drive allele
//		cin>>pa.k;//
//		cin>>pa.epochs;//
//		cin>>pa.compress; //if rand=1, compression level is set at random
		cin>>pa.species;
		cin>>pa.kernel;
		cin>>pa.meanTL;
		cin>>pa.s;
		cin>>pa.EDGEd;
		cin>>pa.singlepop;
		for(int x=0;x<TL;x++)
			{
			cin>>pa.LarvProbs[x]; 
			//cout<<x<<"   "<<pa.LarvProbs[x]<<endl;
			};
		cin>>pa.U;
		cin>>pa.CentRad;


		pa.dispa=0.5;pa.dispb=15;

	//	cout<<pa.U<<endl;
	/*----------------------------------------------------------------------------------------*/
	/*---------------------set mean and variance of local permanent water--------------------------------*/
		if(pa.alpha0>0)
			{
			pa.mu=log(pa.alpha0*pa.alpha0/sqrt(pa.al0var+pa.alpha0*pa.alpha0));
			pa.sig=sqrt(log(1+pa.al0var/(pa.alpha0*pa.alpha0)));	
			};
	/*---------------------------------------------------------------------------------------------------*/

		for(int a=0;a<TL;a++){in.NumJW[a]=10000;};

	/*-------------------------------input rain data file-------------------------------------------------*/
	/*---------------------------------------------------------------------------------------------------*/
	
	/*-------------------------------input mortality-------------------------------------------------*/
	/*---------------------------------------------------------------------------------------------------*/
	/*--------------------------------Set inheritance architecture---------------------------------------*/
		SetFertility();
	/*---------------------------------------------------------------------------------------------------*/

	/*--------------------------------run model NumRuns times--------------------------------------------*/
		RunNReps(ti.NumRuns);
	/*---------------------------------------------------------------------------------------------------*/
	return 0;};
	

	void RunNReps(int N)
	      {

		os3<<"ParList"<<pa.set<<"run"<<pa.index<<".txt";
		ParList.open(os3.str().c_str());
		os1<<"LocalData"<<pa.set<<"run"<<pa.index<<".txt";// make a file for outputing local data
		os2<<"Totals"<<pa.set<<"run"<<pa.index<<".txt";// make a file for outputing global data
		localinfo.open(os1.str().c_str());
		globalinfo.open(os2.str().c_str());
		std::clock_t start;
		double dtime;
		start = std::clock();
		for(int j=0;j<N;j++)
			{
//		pa.xi=rg.Random();
//		SetFertility();
//		ParList<<j+1<<"   "<<pa.xi<<endl;
			initiate();
	for(int pat=0;pat<Site.size();pat+=in.recSitesFreq)ParList<<Site[pat].x<<"   "<<Site[pat].y<<endl;
			RunMaxT();
			};
		dtime = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		os2.str("");
		os1.str("");
		globalinfo.close();
		localinfo.close();
		os3.str("");
		ParList.close();


		return;};

	void RunMaxT(void)
		{
		int TT=ti.start;
		int TTT;
		int swit=0;
		int site;
		int tot,totW,totD,totE,totR;
		int uniquepat,relpat,year2,num_release_clusters,num_driver;
		if(in.NumDriverSites<1)num_release_clusters=int(in.NumDriverSites*to.CentSqVils);
		if(in.NumDriverSites>=1)num_release_clusters=min(to.CentSqVils,int(in.NumDriverSites));
		int relpatches[num_release_clusters];
		int reltimes[num_release_clusters];
		for(int ii=0;ii<num_release_clusters;ii++)relpatches[ii]=500000;
	
		ti.yearnow=0;

		while (TT<ti.maxT+1)
			{
		
			pa.muAF=pa.muA;
			if(TT%365==0)
				{
				for(int i=0;i<Site.size();i++)	Site[i].yearTotF=0; 
				if(num_release_clusters==1){
					for(int i=0;i<Site.size();i++)	if(Site[i].CentSq==1)site=i; 
					relpatches[0]=site;reltimes[0]=in.driver_time;}
				else{
					for(int jj=0;jj<num_release_clusters;jj++) relpatches[jj]=-1;
					for(int jj=0;jj<num_release_clusters;jj++)
						{
						uniquepat=0;
						while(uniquepat==0)
							{
							uniquepat=1;
							relpat=rg.IRandom(0,Site.size()-1);
							for(int ii=0;ii<num_release_clusters;ii++)
								{if(relpat==relpatches[ii] || Site[relpat].CentSq==0)uniquepat=0;};
							relpatches[jj]=relpat;
							};
						if(in.driver_time<0)reltimes[jj]=rg.IRandom(1,364); else reltimes[jj]=in.driver_time;
						};
					};
				};
//		cout<<"B"<<endl;

			if(TT%ti.interval==0)
				{
				totW=0;totD=0;totE=0;totR=0;
				for(int pat=0;pat<Site.size();pat++)
					{
				//	if(Site[pat].CentSq==1)
				//		{
						tot= Site[pat].M[0]+ Site[pat].M[1]+ Site[pat].M[2]+ Site[pat].M[3]+ Site[pat].M[4]+ Site[pat].M[5];
						if(tot==0){totE++;}
						else if(2*Site[pat].M[0]+Site[pat].M[1]+Site[pat].M[3]>tot){totW++;}
						else if(2*Site[pat].M[2]+Site[pat].M[1]+Site[pat].M[5]>tot){totD++;}
						else if(2*Site[pat].M[4]+Site[pat].M[3]+Site[pat].M[5]>tot){totR++;};
				//		};
					};
				cout<<TT;
				for(int i=0;i<NumGen;i++)cout<<"     "<<to.Em[i];
			//	for(int i=0;i<NumGen;i++)cout<<"     "<<to.M[i];
				cout<<"   "<<totW<<"   "<<totD<<"    "<<totR<<"    "<<totE<<"     "<<to.CentSqVils<<endl;
				//globalinfo<<TT<<"   "<<to.Fww<<"     "<<to.Fwd<<"     "<<to.Fdd<<"     "<<to.Fwr<<"     "<<to.Frr<<"     "<<to.Fdr<<endl;
			for(int i=0;i<NumGen;i++)to.CentF[i]=0;
	for(int pat=0;pat<Site.size();pat++)
		{
		if(Site[pat].CentSq==1)
			{
			for(int i=0;i<NumGen;i++)to.CentF[i]+=accumulate(Site[pat].F[i],Site[pat].F[i]+NumGen,0);
			};
		};
				globalinfo<<TT;
				for(int i=0;i<NumGen;i++)globalinfo<<"     "<<to.CentF[i];
				globalinfo<<endl;
				};
			TT++;

				for(int i=0;i<NumGen;i++)to.Em[i]=0;
			//	cout<<TT<<" C1   "<<num_release_clusters<<"     "<<reltimes[0]<<"    "<<relpatches[0]<<"    "<<to.CentSqVils<<"    "<<in.NumDriverSites<<endl;
			if(TT>=abs(in.driver_start) && TT<in.driver_end)
			{
		//		cout<<TT<<" A   "<<num_release_clusters<<"     "<<reltimes[0]<<"    "<<relpatches[0]<<endl;
				for(int jj=0;jj<num_release_clusters;jj++)
					{
		//		cout<<TT<<"  B  "<<num_release_clusters<<endl;
					if(TT%365==reltimes[jj])
						{
			//	cout<<TT<<" C2   "<<jj<<"    "<<num_release_clusters<<"     "<<reltimes[jj]<<"    "<<relpatches[jj]<<endl;
						PutDriverSites(relpatches[jj]);
						};
					};
			};
			/*TTT=TT-abs(in.driver_time)+365;
			if(TTT%ti.rec==292){record(TT);};*/
			if(TT>ti.recstart && TT<=ti.recend && TT%ti.recfreq==0)record(TT);
			//if(TT==ti.maxT)record(TT);
			to.distD=0;
			OneStep(TT);
			};
        return;};


					
	void initiate(void){
		for(int i=0;i<NumGen;i++)
		{
		to.Em[i]=0; 
		to.J[i]=0;
		to.M[i]=0;
		to.V[i]=0;
		to.F[i]=0;
		to.CentF[i]=0;
		};
		to.JTot=0;to.MTot=0;to.VTot=0;to.FTot=0;
		to.distW=0;to.distD=0;to.distR=0;
		to.CentSqVils=0;
		to.CentSqHum=0;
		Site.clear();
		for(int xx=0;xx<nx;xx++) { for(int yy=0;yy<ny;yy++) { SetsPerCell[xx][yy].clear(); }; };
		Patch pp;
		//cout<<"init   1"<<endl;
		/*----------------------------------------------------------------------------------*/
		/*------------------------------input the settlement data---------------------------*/
		for(int ii=0;ii<pa.NumPat;ii++)
		{       
			pp.x=rg.Random()*pa.U;
			pp.y=rg.Random()*pa.U;
			pp.WaterTemp=0;
			pp.WaterPerm=0;
			pp.sqx=0;
			pp.sqy=0;
			pp.CentSq=0;
			if(dist(pp.x,pp.y,pa.U/2.0,pa.U/2.0)<pa.CentRad)pp.CentSq=1;

			pp.Nhum=1;
			pp.gam=0;
			pp.arab=0;
			pp.fun=0;
			pp.area=1;
			pp.connecIND.clear();
			pp.connecW.clear();
			pp.TotW=0;
			//cout<<pp.x<<"     "<<pp.y<<endl;
			for(int i=0;i<NumGen;i++)
			{
			for(int a=0;a<TL;a++)
				{
				pp.J[i][a]=0; 
				};
			pp.M[i]=0;
			pp.V[i]=0;
			for(int j=0;j<NumGen;j++)
			{
			pp.F[i][j]=0;
			pp.AesF[i][j]=0;
			};
			};
			pp.comp=0;pp.mate_rate=0;pp.JTot=0;pp.MTot=0;
			pp.distOrigin=0;
			pp.dens=0.0000817992*pp.Nhum/pp.area;//scaler converts degrees to km^2
			//cout<<pp.Nhum<<"   "<<pp.area<<"    "<<pp.dens<<endl;
			if(pa.alpha0>0)pp.alpha0=exp(random_normal(pa.mu,pa.sig));else pp.alpha0=0;
			if(pp.CentSq==1){to.CentSqVils++;to.CentSqHum+=pp.Nhum;};
			Site.push_back(pp);
			//SetsPerCell[pp.sqx][pp.sqy].push_back(Site.size()-1);
		};
		/*----------------------------------------------------------------------------------*/

		//kMeansClustering();

		//cout<<"numdriver   "<<in.NumDriverClus<<"    "<<in.NumDriver<<endl;
//		kmeans<<pa.compress<<"    "<<pa.k<<"     "<<Site.size()<<"    "<<Site.size()<<endl;
//		cout<<"number of sites:   "<<Site.size()<<"  num of sites in centre  "<<to.CentSqVils<<endl;
	//	exit(1);	
	/*-------------------------------numbers to use when initiating populations--------------------------*/
		int inV=10000;
		int inM=50000;
		int inF=40000;
		in.NumAdultsWV=inV; in.NumAdultsWM=inM; in.NumAdultsWF=inF;
		in.NumAdultsRV=inV; in.NumAdultsRM=inM; in.NumAdultsRF=inF;
		in.NumAdultsDV=inV; in.NumAdultsDM=inM; in.NumAdultsDF=inF;
	/*---------------------------------------------------------------------------------------------------*/
		//cout<<Site.size()<<endl;
		for(int i=0;i<Site.size();i++){if(Site[i].CentSq==1)SitesPopulate(i);};

		UpdateConnec();
//cout<<"finished initialising"<<endl;
	return;};


	void UpdateConnec()
	{
	double dd,ww,cor;
	if(pa.kernel==-1)
	{
	cor=((pa.dispb-2)*(pa.dispb-1)*(pa.dispa*pa.dispa-(pa.dispa+pa.LD)*std::pow((pa.dispa+pa.LD)/pa.dispa,-1*pa.dispb)*(pa.dispa+(pa.dispb-1)*pa.LD)))/(pa.dispa*pa.dispa*(2-3*pa.dispb+pa.dispb*pa.dispb));
	//cor=(std::pow((pa.disp_a + pa.LD)/ pa.disp_a),-1*pa.disp_b)* (-pa.disp_a*pa.disp_b*pa.LD - (-1 + pa.disp_b)*pa.LD^2 + pa.disp_a^2 (-1 + ((pa.disp_a + pa.LD)/pa.disp_a)^pa.disp_b)))/pa.disp_a^2;
	for(int index=0;index<Site.size();index++)
	{
		Site[index].connecIND.clear();
		Site[index].connecW.clear();
		Site[index].TotW=0;
		for(int ii=0;ii<Site.size();ii++)
			{
				dd=dist(Site[index].x,Site[index].y,Site[ii].x,Site[ii].y);
				if(dd<pa.LD)
				{
					Site[index].connecIND.push_back(ii); 
					//ww=((((pa.disp_a - 1)*(pa.disp_b - 2)/(2*Pi*pa.disp_a^2)) (1 + dd/pa.disp_a)^-pa.disp_b))/cor;
					ww=((((pa.dispb - 1)*(pa.dispb - 2)/(2*PI*std::pow(pa.dispa,2)))*std::pow(1 + dd/pa.dispa,-1*pa.dispb)))/cor;
					Site[index].connecW.push_back(ww); 
					Site[index].TotW+=ww;
				};
			};
		//	cout<<Site[index].connecIND.size()<<"   "<<Site[index].TotW<<endl;

	};
	};
	if(pa.kernel==0)
	{
	for(int index=0;index<Site.size();index++)
	{
		Site[index].connecIND.clear();
		Site[index].connecW.clear();
		Site[index].TotW=0;
		for(int ii=0;ii<Site.size();ii++)
			{
			if(ii!=index)
				{
				//dd=distB(1.0,Site[index].x,Site[index].y,Site[ii].x,Site[ii].y);
				dd=distB(pa.U,Site[index].x,Site[index].y,Site[ii].x,Site[ii].y);
			//	dd=dist(Site[index].x,Site[index].y,Site[ii].x,Site[ii].y);
				if(dd<pa.LD)
				{
					Site[index].connecIND.push_back(ii); 
					ww=(1-dd/pa.LD)*(3/(PI*pa.LD*pa.LD))/(1.0*pa.NumPat);
					Site[index].connecW.push_back(ww); 
					Site[index].TotW+=ww;
				};
				};
			};
		//	cout<<Site[index].connecIND.size()<<"   "<<Site[index].TotW<<endl;

	};
	};


	return;};


	void PutDriverSites(int pat){
		long long int num_driver=(long long int) in.NumDriver*(1.0*Site.size())/(1.0*Site.size());
		long long int driver=random_poisson(1.0*num_driver);
		Site[pat].M[1]+=driver;
		Site[pat].MTot+=driver;
		to.M[1]+=driver;
		to.MTot+=driver;
		UpdateMate();
		return;};


	void SitesPopulate(int pat){
		//	for(int a=0;a<TL;a++){Site[pat].J[0][a]+=in.NumJW[a];to.J[0]+=in.NumJW[a];to.JTot+=in.NumJW[a];Site[pat].JTot+=in.NumJW[a];};
			Site[pat].M[0]=int(in.NumAdultsWM);to.M[0]+=int(in.NumAdultsWM);to.MTot+=int(in.NumAdultsWM);Site[pat].MTot+=int(in.NumAdultsWM);
			Site[pat].V[0]=int(in.NumAdultsWV);to.V[0]+=int(in.NumAdultsWV);to.VTot+=int(in.NumAdultsWV);
			Site[pat].F[0][0]=int(in.NumAdultsWF);
			to.F[0]+=int(in.NumAdultsWF);
			to.FTot+=int(in.NumAdultsWF);
			if(Site[pat].distOrigin>to.distW)to.distW=Site[pat].distOrigin;
			Site[pat].comp=std::pow((double)pa.alpha0/(pa.alpha0+Site[pat].JTot),1.0/pa.meanTL);
			Site[pat].mate_rate=Site[pat].MTot/(pa.beta+Site[pat].MTot);

//cout<<"SP  "<<pat<<"   "<<Site[pat].CentSq<<"    "<<Site[pat].F[0][0]<<"   "<<Site[pat].J[0][0]<<"   "<<Site[pat].comp<<"   "<<to.FTot<<endl;

	return;};
	void OneStep(int day){
		JuvGetOlder();
		AdultsDie();
		VirginsMate();
		AdultsMove();
		LayEggs();
		JuvEmerge();
		if(day%365 > pa.t_hide1 && day%365<=pa.t_hide2 && pa.psi>0.00001)Hide();
		if(day%365 > pa.t_wake1 && day%365<=pa.t_wake2 &&pa.psi>0.00001)Wake(day);
		UpdateComp(day);
		UpdateMate();
	return;};
	

	void JuvEmerge(void){
		int surv,survM;
		double pcomp;
		for(int pat=0;pat<Site.size();pat++)
		{
		pcomp=(1-pa.muJ)*Site[pat].comp;
		for(int i=0;i<NumGen;i++)
		{
			surv=random_binomial(Site[pat].J[i][0],pcomp); 
			if(Site[pat].CentSq==1)to.Em[i]+=surv;
			if(surv>0)
			{		
				survM=random_binomial(surv,0.5); Site[pat].MTot+=survM;
				Site[pat].M[i]+=survM; 
				to.M[i]+=survM;
				to.MTot+=survM;
				Site[pat].V[i]+=surv-survM;
				to.V[i]+=surv-survM;
				to.VTot+=surv-survM;
			if((i==1||i==2||i==5) && Site[pat].distOrigin>to.distD){ to.distD=Site[pat].distOrigin;};
			};
		};
		};
		return;};

	void JuvGetOlder(void){
		long int jtSitesch,jtAllAll;
		long int jtAll[NumGen];
		for(int i=0;i<NumGen;i++) { jtAll[i]=0; };
		double pcomp;
		for(int pat=0;pat<Site.size();pat++)
		{
		pcomp=(1-pa.muJ)*Site[pat].comp;
			jtSitesch=0;
			for(int age=0;age<TL-1;age++)
			{
				for(int i=0;i<NumGen;i++)
					{		
					Site[pat].J[i][age]=random_binomial(Site[pat].J[i][age+1],pcomp);
					jtSitesch+=Site[pat].J[i][age];
					jtAll[i]+=Site[pat].J[i][age];
					};
			};
			
			Site[pat].JTot=jtSitesch;
			jtAllAll+=jtSitesch;
		for(int i=0;i<NumGen;i++)Site[pat].J[i][TL-1]=0; 
		};
		for(int i=0;i<NumGen;i++)to.J[i]=jtAll[i];
		to.JTot=jtAllAll;
	return;};

	void VirginsMate(){
		int v[NumGen];
		int f[NumGen];
		int* mates;
		int num;
		for(int pat=0;pat<Site.size();pat++)
				{
			for(int i=0;i<NumGen;i++)
					{
						v[i]=random_binomial(Site[pat].V[i],Site[pat].mate_rate);
					if(v[i]>0)
					{
					mates=random_multinom(v[i],Site[pat].M);
					for(int j=0;j<NumGen;j++)Site[pat].F[i][j]+=*(mates+j);
					Site[pat].yearTotF+=v[i];
					to.F[i]+=v[i];
					delete[] mates;
					Site[pat].V[i]-=v[i]; to.V[i]-=v[i]; to.VTot-=v[i]; to.FTot+=v[i];
					};
					};
				};
	return;};






	void AdultsMove(void){
		int pat,newpat;
		int* aww;
		int howmany;
		double pdisperse;
		if(Site.size()>1)
		{
		for(pat=0;pat<Site.size();pat++)
			{
		if(pa.kernel==-1) { pdisperse=pa.d*(1-exp(-1*Site[pat].TotW)); };
		//if(pa.kernel==0) { pdisperse=pa.d*(1-exp(-1*Site[pat].TotW)); };
		if(pa.kernel==0) { pdisperse=pa.d; };
//cout<<"pdisp  "<<pdisperse<<"     "<<Site[pat].area<<endl;
				for(int i=0;i<NumGen;i++)
				{
					Site[pat].MoveM[i]=random_binomial(Site[pat].M[i],pdisperse);
					Site[pat].M[i]-=Site[pat].MoveM[i];
					Site[pat].MTot-=Site[pat].MoveM[i];
					for(int j=0;j<NumGen;j++)
						{
						Site[pat].MoveF[i][j]=random_binomial(Site[pat].F[i][j],pdisperse);
						Site[pat].F[i][j]-=Site[pat].MoveF[i][j]; 
						};
				};
			};
		for(pat=0;pat<Site.size();pat++)
				{
				howmany=Site[pat].connecIND.size();
				for(int i=0;i<NumGen;i++)
				{
					if(Site[pat].MoveM[i]>0)
						{
						aww=random_multinom_var(Site[pat].MoveM[i],howmany,&Site[pat].connecW[0],Site[pat].TotW);
						for(newpat=howmany-1;newpat>=0;newpat--)
							{
								Site[Site[pat].connecIND[newpat]].M[i]+=*(aww+newpat);
								Site[Site[pat].connecIND[newpat]].MTot+=*(aww+newpat);
							};
						delete[] aww;
						};
					for(int j=0;j<NumGen;j++)
						{
						if(Site[pat].MoveF[i][j]>0)
							{
							aww=random_multinom_var(Site[pat].MoveF[i][j],howmany,&Site[pat].connecW[0],Site[pat].TotW);
							for(newpat=howmany-1;newpat>=0;newpat--) Site[Site[pat].connecIND[newpat]].F[i][j]+=*(aww+newpat);
							delete[] aww;
							};
						};
				};
				};
		};
	return;};



	void Hide(){
		long long int num;
		for(int pat=0;pat<Site.size();pat++)
			{
			for(int i=0;i<NumGen;i++)
				{
				for(int j=0;j<NumGen;j++)
					{
					num=random_binomial(Site[pat].F[i][j],pa.psi);Site[pat].F[i][j]-=num;
					to.F[i]-=num;
					to.FTot-=num;
					Site[pat].AesF[i][j]+=random_binomial(num,1-pa.muAES);	
					};
				};
			};
		return;};

	void Wake(int day){
		long long int num;
		double prob=1.0/(1.0+pa.t_wake2-(day%365));
		for(int pat=0;pat<Site.size();pat++)
				{
				for(int i=0;i<NumGen;i++)
					{
					for(int j=0;j<NumGen;j++)
						{
						num=random_binomial(Site[pat].AesF[i][j],prob);
						Site[pat].F[i][j]+=num;
						to.F[i]+=num;
						to.FTot+=num;
						Site[pat].AesF[i][j]-=num;	
						};
					};
				};
		return;};




	void LayEggs (void){
		double num;
		int* aww;
		for(int pat=0;pat<Site.size();pat++)
				{
	//fraction of 0:ww,1:wd,2:dd,3:wr,4:rr,5:dr
				for(int i=0;i<NumGen;i++)
					{
					for(int j=0;j<NumGen;j++)
						{
						for(int k=0;k<NumGen;k++)
							{
							num=Site[pat].F[i][j]*pa.f[i][j][k];
					//		cout<<num<<"   ";
							num=random_poisson(pa.theta*num);
							aww=random_multinom_var(num,TL,&pa.LarvProbs[0],1.0);
							for(int a=0;a<TL;a++) Site[pat].J[k][a]+=*(aww+a);
							delete[] aww;
							to.J[k]+=num;
							to.JTot+=num;
							};
						};
					};

				};
		return;};


	void AdultsDie(void){
		long long int num;
		for(int pat=0;pat<Site.size();pat++)
				{
				for(int i=0;i<NumGen;i++)
					{
					num=random_binomial(Site[pat].M[i],pa.muA);
					Site[pat].M[i]-=num;
					Site[pat].MTot-=num;
					to.M[i]-=num;
					to.MTot-=num;	
					num=random_binomial(Site[pat].V[i],pa.muA);
					Site[pat].V[i]-=num;
					to.V[i]-=num;
					to.VTot-=num;	

					for(int j=0;j<NumGen;j++)
						{
						num=random_binomial(Site[pat].F[i][j],pa.muAF);
						Site[pat].F[i][j]-=num;
						to.F[i]-=num;
						to.FTot-=num;	
						};
					};
				};
		return;};

void UpdateComp(int day){
	int week=(int)(day % 365)/7;
	if(week==52)week=51;
	double al;
	for(int pat=0;pat<Site.size();pat++)
		{
		al=0;
	/*-----unnormalised human pop scaling---------*/
			al=Site[pat].alpha0;
	/*-------------------------------------------*/
			Site[pat].comp=std::pow(al/(al+Site[pat].JTot+0.0001),1/pa.meanTL);
//		if(pat==0)cout<<Site[0].comp<<"    "<<al<<endl;
		};
		return;};

	void UpdateMate(void){
		for(int pat=0;pat<Site.size();pat++)
				{
				Site[pat].mate_rate=Site[pat].MTot/(pa.beta+Site[pat].MTot);
				};
		return;};

	void SetFertility()
{
//fraction of 0:ww,1:wd,2:dd,3:wr,4:rr,5:dr
	double Fwwww[6]={1,0,0,0,0,0};
	double Fwwwd[6]={(1-pa.e-pa.gamma)*0.5,(1+pa.e)*0.5,0,pa.gamma*0.5,0,0};
	double Fwwdd[6]={0,1,0,0,0,0};
	double Fwwwr[6]={0.5,0,0,0.5,0,0};
	double Fwwrr[6]={0,0,0,1,0,0};
	double Fwwdr[6]={0,0.5,0,0.5,0,0};
	double Fwdww[6]={(1-pa.xi)*(1-pa.e-pa.gamma)*0.5,(1-pa.xi)*(1+pa.e)*0.5,0,(1-pa.xi)*pa.gamma*0.5,0,0};
	double Fwdwd[6]={(1-pa.xi)*(1-pa.e-pa.gamma)*(1-pa.e-pa.gamma)*0.25,(1-pa.xi)*(1-pa.e-pa.gamma)*(1+pa.e)*0.5,(1-pa.xi)*(1+pa.e)*(1+pa.e)*0.25,(1-pa.xi)*(1-pa.e-pa.gamma)*pa.gamma*0.5,(1-pa.xi)*pa.gamma*pa.gamma*0.25,(1-pa.xi)*(1+pa.e)*pa.gamma*0.5};

	double Fwddd[6]={0,(1-pa.xi)*(1-pa.e-pa.gamma)*0.5,(1-pa.xi)*(1+pa.e)*0.5,0,0,(1-pa.xi)*pa.gamma*0.5};

	double Fwdwr[6]={(1-pa.xi)*(1-pa.e-pa.gamma)*0.25,(1-pa.xi)*(1+pa.e)*0.25,0,(1-pa.xi)*((1-pa.e-pa.gamma)*0.25+(pa.gamma*0.25)),(1-pa.xi)*pa.gamma*0.25,(1-pa.xi)*(1+pa.e)*0.25};

	double Fwdrr[6]={0,0,0,(1-pa.xi)*(1-pa.e-pa.gamma)*0.5,(1-pa.xi)*pa.gamma*0.5,(1-pa.xi)*(1+pa.e)*0.5};
	double Fwddr[6]={0,(1-pa.xi)*(1-pa.e-pa.gamma)*0.25,(1-pa.xi)*(1+pa.e)*0.25,(1-pa.xi)*(1-pa.e-pa.gamma)*0.25,(1-pa.xi)*pa.gamma*0.25,(1-pa.xi)*((1+pa.e)*0.25+pa.gamma*0.25)};
	double Fddww[6]={0,0,0,0,0,0};
	double Fddwd[6]={0,0,0,0,0,0};
	double Fdddd[6]={0,0,0,0,0,0};
	double Fddwr[6]={0,0,0,0,0,0};
	double Fddrr[6]={0,0,0,0,0,0};
	double Fdddr[6]={0,0,0,0,0,0};
	double Fwrww[6]={0.5,0,0,0.5,0};
	double Fwrwd[6]={(1-pa.e-pa.gamma)*0.25,(1+pa.e)*0.25,0,(pa.gamma*0.25+(1-pa.e-pa.gamma)*0.25),pa.gamma*0.25,(1+pa.e)*0.25};
	double Fwrdd[6]={0,0.5,0,0,0,0.5};
	double Fwrwr[6]={0.25,0,0,0.5,0.25,0};
	double Fwrrr[6]={0,0,0,0.5,0.5,0};
	double Fwrdr[6]={0,0.25,0,0.25,0.25,0.25};
	double Frrww[6]={0,0,0,0,0,0};
	double Frrwd[6]={0,0,0,0,0,0};
	double Frrdd[6]={0,0,0,0,0,0};
	double Frrwr[6]={0,0,0,0,0,0};
	double Frrrr[6]={0,0,0,0,0,0};
	double Frrdr[6]={0,0,0,0,0,0};

	double Fdrww[6]={0,0,0,0,0,0};
	double Fdrwd[6]={0,0,0,0,0,0};
	double Fdrdd[6]={0,0,0,0,0,0};
	double Fdrwr[6]={0,0,0,0,0,0};
	double Fdrrr[6]={0,0,0,0,0,0};
	double Fdrdr[6]={0,0,0,0,0,0};




	for(int k=0;k<6;k++)
	{
		for(int i=0;i<6;i++)
		{
			for(int j=0;j<6;j++)
			{
				if(i==0)
				{
					if(j==0) pa.f[i][j][k]=Fwwww[k];
					if(j==1) pa.f[i][j][k]=Fwwwd[k];
					if(j==2) pa.f[i][j][k]=Fwwdd[k];
					if(j==3) pa.f[i][j][k]=Fwwwr[k];
					if(j==4) pa.f[i][j][k]=Fwwrr[k];
					if(j==5) pa.f[i][j][k]=Fwwdr[k];
				};
				if(i==1)
				{
					if(j==0) pa.f[i][j][k]=Fwdww[k];
					if(j==1) pa.f[i][j][k]=Fwdwd[k];
					if(j==2) pa.f[i][j][k]=Fwddd[k];
					if(j==3) pa.f[i][j][k]=Fwdwr[k];
					if(j==4) pa.f[i][j][k]=Fwdrr[k];
					if(j==5) pa.f[i][j][k]=Fwddr[k];
				};
				if(i==2)
				{
					if(j==0) pa.f[i][j][k]=Fddww[k];
					if(j==1) pa.f[i][j][k]=Fddwd[k];
					if(j==2) pa.f[i][j][k]=Fdddd[k];
					if(j==3) pa.f[i][j][k]=Fddwr[k];
					if(j==4) pa.f[i][j][k]=Fddrr[k];
					if(j==5) pa.f[i][j][k]=Fdddr[k];
				};
				if(i==3)
				{
					if(j==0) pa.f[i][j][k]=Fwrww[k];
					if(j==1) pa.f[i][j][k]=Fwrwd[k];
					if(j==2) pa.f[i][j][k]=Fwrdd[k];
					if(j==3) pa.f[i][j][k]=Fwrwr[k];
					if(j==4) pa.f[i][j][k]=Fwrrr[k];
					if(j==5) pa.f[i][j][k]=Fwrdr[k];
				};
				if(i==4)
				{
					if(j==0) pa.f[i][j][k]=Frrww[k];
					if(j==1) pa.f[i][j][k]=Frrwd[k];
					if(j==2) pa.f[i][j][k]=Frrdd[k];
					if(j==3) pa.f[i][j][k]=Frrwr[k];
					if(j==4) pa.f[i][j][k]=Frrrr[k];
					if(j==5) pa.f[i][j][k]=Frrdr[k];
				};
				if(i==5)
				{
					if(j==0) pa.f[i][j][k]=Fdrww[k];
					if(j==1) pa.f[i][j][k]=Fdrwd[k];
					if(j==2) pa.f[i][j][k]=Fdrdd[k];
					if(j==3) pa.f[i][j][k]=Fdrwr[k];
					if(j==4) pa.f[i][j][k]=Fdrrr[k];
					if(j==5) pa.f[i][j][k]=Fdrdr[k];
				};
			};
		};
	};	

	//for(int k=0;k<6;k++) { for(int i=0;i<6;i++) { for(int j=0;j<6;j++) { cout<<i<<"   "<<j<<"    "<<k<<"   "<<pa.f[i][j][k]<<"    "; };cout<<endl; }; };
	
	

	return;};

	long long int longmin(long long int a,long long int b)
			{
			long long int c=b;
			if(a<b)c=a;
			return c;
			};

	long long int longmax(long long int a,long long int b)
			{
			long long int c=b;
			if(a>b)c=a;
			return c;
			};

	long long int random_binomial(long long int N,double p){
		long long int ran;
		if(N==0){ran=0;}
		else if(p>0.999999){ran=N;}
		else if(p<0.000001){ran=0;}
		else if(N*p>10 && N*(1-p)>10) {ran=longmax(0,longmin(N,(long long int)random_normal(N*p,sqrt(N*p*(1-p)))));}
		else if((N>20 && p<0.05) || (N>100 && N*p<10)){ran=random_poisson(N*p);}
		else if((N>20 && p>0.95) || (N>100 && N*(1-p)<10)){ran=N-random_poisson(N*(1-p));}
		else { variate_generator<mt19937, binomial_distribution<> > b(mt19937(int(rg.Random()*time(NULL))), binomial_distribution<>(N,p)); ran=b();};
		return ran;};


       double random_normal(double mu, double sigma)
{
	const double epsilon = std::numeric_limits<double>::min();

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;

	double u1, u2;
	do
	 {
	   u1 = rg.Random();
	   u2 = rg.Random();
	 }
	while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(TWOPI * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(TWOPI* u2);
	return z0 * sigma + mu;
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

	long long int random_poisson(double landa)
	{
		long long int k;
		if(landa<1e-5){k=0;}
		else if(landa>30) {k=longmax(0,(long long int)random_normal(landa,sqrt(landa)));}
		else
			{
			double p=exp(-landa);
			double g=p;
			double u=rg.Random()*0.999999999;
			k=0;
			while (u>g)
			    {
				p*=(landa/(double)(++k));
				g+=p;
			};
			};
	return k;
	};

	int* random_multinom(int N,long long int probs[6])
{
	int *ran=new int[6];
	double sum=double(probs[0]+probs[1]+probs[2]+probs[3]+probs[4]+probs[5]);
	ran[0]=random_binomial(N,probs[0]/sum);
	ran[1]=random_binomial(N-ran[0],probs[1]/(sum-probs[0]));
	ran[2]=random_binomial(N-ran[0]-ran[1],probs[2]/(sum-probs[0]-probs[1]));
	ran[3]=random_binomial(N-ran[0]-ran[1]-ran[2],probs[3]/(sum-probs[0]-probs[1]-probs[2]));
	ran[4]=random_binomial(N-ran[0]-ran[1]-ran[2]-ran[3],probs[4]/(sum-probs[0]-probs[1]-probs[2]-probs[3]));
	ran[5]=N-ran[0]-ran[1]-ran[2]-ran[3]-ran[4];
	return ran;};

int* random_multinom_var(int N,int howmany,double *relprobs,double tot)
{
	int *ran=new int[howmany];
	double sum=tot;
	int Nused=N;
	for(int pat=0;pat<howmany;pat++)
	{
		if(N>0)
		{
		ran[pat]=random_binomial(Nused,*(relprobs+pat)/sum);
		sum-=*(relprobs+pat);
		Nused-=ran[pat];
		}
		else ran[pat]=0;
	};
	return ran;};

int* random_multinomEqualProb(int N,int howmany)
{
	double prob=1/(1.0*howmany);
	int *ran=new int[howmany];
	double sum=1;
	int Nused=N;
	for(int pat=0;pat<howmany;pat++)
	{
		if(N>0)
		{
		ran[pat]=random_binomial(Nused,prob/sum);
		sum-=prob;
		Nused-=ran[pat];
		}
		else ran[pat]=0;
	};
	return ran;};


void record(int index)
	{	
	char state;
	//int tot; 
	for(int pat=0;pat<Site.size();pat+=in.recSitesFreq)
				{
				//localinfo<<index<<"    "<<Site[pat].x<<"   "<<Site[pat].y;
				for(int i=0;i<NumGen;i++) localinfo<<Site[pat].M[i]<<"    "; localinfo<<Site[pat].M[5]<<endl;
				/*tot= Site[pat].M[0]+ Site[pat].M[1]+ Site[pat].M[2]+ Site[pat].M[3]+ Site[pat].M[4]+ Site[pat].M[5];
				if(tot==0){state='E';	}
				else if(2*Site[pat].M[0]>tot){state='W';}
				else if(2*Site[pat].M[2]+Site[pat].M[1]>tot)state='D';
				localinfo<<state<<endl;*/	
			
	//if(Site[pat].CentSq==1)localinfo<<Site[pat].yearTotF<<endl;
				};

		
	return;};




//Below is code for the random number generator---------------------------------------------------------------------------------------------------------------------------------



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
uint32 CRandomMersenne::BRandom() {// Generate 32 random bits

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
   union {float f; uint32 i[2];} convert; //Union allows one portion of the memory to be accessed as different data types.
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


