#include<mpi.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include</home/joyjit/mt19937ar.h>
#include</home/joyjit/mt19937ar.c>

/* This prog simulates the long rods problem with no intersection
 * allowed in two dimensions. The rods are monodiepersed.
 * It is a montecarlo simulation with heat
 * bath dynamics where all the horizontal k-mers are evaporated and
 * redeposited, then vertical k-mers and so on.  The initial
 * conditions can be of two kinds: (1) all empty and (2) half vertical
 * half horizontal. This is controlled using INITIAL_FLAG. The
 * mean density, order parameter and second and fourth moment are
 * measured. Horizontal rod is marked with 1 on lat while vertical rod
 * is marked with 2. */

#define K 7  /* length of K-mer*/
#define L 1680	/*Number of sites in one direction, should be even*/
#define N (L*L) /* The total number of sites */
#define DIR 2 /*Number of positive directions */
#define NN 2 /*Number of nearest neighbours in forward direction */
#define T 5000000	/*# of MC steps for equilibration*/
#define BLOCKS 10 /*Number of blocks in which average is taken*/
#define AVE 20000000 /* # of readings in each block */
#define PSTART 0.7 /*starting probability*/
#define PDIFF 0.1 /*increment in probability*/
#define PEND (PSTART+0.5*PDIFF) /*ending probability */
#define GAP2 25000 /*GAP after which output is written during equilibration */
#define INITIAL_FLAG 0 /*0: all empty 1: all filled */
#define BINS 101 /*Keep it odd*/
#define BINSIZE (2.0/(1.0*BINS))
#define MASTER 0               /* taskid of first task */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW(id,p,n))

int lasth,lastv,lasth0,lastv0;
int starth,finalh,endh,startv,finalv,endv,taskid,numworkers,numtasks,rows,columns;
double count[BINS],periodic[K],acceptance[L+1],mu;
double meanrho[BLOCKS],meanabs[BLOCKS],rho2[BLOCKS],fluc[BLOCKS];
double meanm1[BLOCKS],meanm2[BLOCKS],meanm4[BLOCKS];
double Q_1[BLOCKS],Q_2[BLOCKS],Q_4[BLOCKS],Rho_2[BLOCKS],Rho[BLOCKS],modQ[BLOCKS];
double mass[BINS],factor[BINS],bincount[BINS];
char outfile1[100],outfile2[100],outfile3[100];
char readfile1[100],readfile2[100];
long seedval;
double elapsed_time;
int *lat;
int *ln;
int *rn;
int *bn;
int *tn;
int *horizontal;
int *vertical;
int *displs;
int *scounts;
MPI_Datatype columntype;
MPI_Status status;
void print_config()
{
	/* Prints out the configuration */
	int x,y;

	for(y=L-1;y>=0;y--)
	{
		for(x=0;x<L;x++)
			printf("%2d",lat[find_site(x,y)]);
		printf("\n");
	}
	printf("***************************\n");
	sleep(2);
}

int find_site(int x, int y)
{
	/* Given x, y it gives site index */
	x=(x+L) % L;
	y=(y+L) % L;
	return x+y*L;
}

/*void take_input()
{
	seedval=(taskid+2)*825-4*taskid;
	init_genrand(seedval);
}*/

long seedgen(void)
{
	long s, seed, pid, seconds;
	pid = getpid();
	s = time ( &seconds );
	seed = abs(((s*181)*((pid-83)*359))%103929);
	return seed;
}

void take_input()
{
	seedval=seedgen();
	init_genrand(seedval);
}

void neighbour()
{
	int i,j,x,y;
	for(y=0;y<L;y++)
	{
		for(x=0;x<L;x++)
		{
			i=find_site(x,y);
			j=find_site(x-1,y); ln[i]=j;
			j=find_site(x+1,y); rn[i]=j;
			j=find_site(x,y-1); bn[i]=j;
			j=find_site(x,y+1); tn[i]=j;
		}
	}
}

void initialize()
{
	/* initializes nbr list and output file names */
	int i,j,x,y;
	double tmp,tmp1;

	/*for(y=0;y<L;y++)
	{
		for(x=0;x<L;x++)
		{
			i=find_site(x,y);
			j=find_site(x-1,y); ln[i]=j;
			j=find_site(x+1,y); rn[i]=j;
			j=find_site(x,y-1); bn[i]=j;
			j=find_site(x,y+1); tn[i]=j;
		}
	}*/
	if(INITIAL_FLAG==0)
		sprintf(outfile1,"emptyK%dL%dM%5.4lf",K,L,mu);
	if(INITIAL_FLAG==1)
		sprintf(outfile1,"filledK%dL%dM%5.4lf",K,L,mu);
	/*initializing bin parameters*/

	for(j=0;j<BINS;j++)
	{
		bincount[j]=0;
		mass[j]=0;
	}
	for(j=-L*L;j<=L*L;j=j+K)
	{
		tmp=(1.0*j)/(1.0*N);
		i=floor((tmp+1.0)/BINSIZE);
		if(j==-N)
			i=0;
		if(j==N)
			i=BINS-1;
		bincount[i]++;
		mass[i]=mass[i]+tmp;
	}
	tmp1=(K*1.0)/(1.0*N);
	for(j=0;j<BINS;j++)
	{
		mass[j]=mass[j]/bincount[j];
		factor[j]=tmp1*bincount[j];
	}
}

void lat_initialization()
{
	int i;
	lastv=0;lasth=0; 
	if(INITIAL_FLAG==0)
	{
		for(i=0;i<N;i++)
			lat[i]=0;
	}
	
	else if(INITIAL_FLAG==1)
	{
		for(i=0;i<N/2;i++)
		{
			lat[i]=1;
			horizontal[lasth]=i;
			lasth++;
		}
		for(i=N/2;i<N;i++)
		{
			lat[i]=2;
			vertical[lastv]=i;
			lastv++;
		}
	}
}

void lat_init()
{
	/* Initializes quantities that have to initialized for every value
	 * of probability p */
	int i;
	double x;
	FILE *fp;
	/*lastv=0;lasth=0;*/ 
	for(i=0;i<BLOCKS;i++)
	{
		meanrho[i]=0;meanabs[i]=0;meanm1[i]=0; meanm2[i]=0; meanm4[i]=0,rho2[i]=0,fluc[i]=0;
		Rho[i]=0;Rho_2[i]=0;modQ[i]=0;Q_1[i]=0;Q_2[i]=0;Q_4[i]=0;
	}
	if((L % K !=0)|| (L % 2) !=0)
	{
		printf("ERROR IN DIVISIBILITY\n");
		exit(0);
	}
	for(i=0;i<BINS;i++)
		count[i]=0;
	if(INITIAL_FLAG==0)
	{
		/*for(i=0;i<N;i++)
			lat[i]=0;*/
		sprintf(outfile2,"emptyK%dL%dM%5.4lf_t",K,L,mu);
		sprintf(outfile3,"emptyK%dL%dM%5.4lf_p",K,L,mu);
	}
	if(INITIAL_FLAG==1)
	{
		/*for(i=0;i<N/2;i++)
		{
			lat[i]=1;
			horizontal[lasth]=i;
			lasth++;
		}
		for(i=N/2;i<N;i++)
		{
			lat[i]=2;
			vertical[lastv]=i;
			lastv++;
		}*/
		sprintf(outfile2,"filledK%dL%dM%5.4lf_t",K,L,mu);
		sprintf(outfile3,"filledK%dL%dM%5.4lf_p",K,L,mu);
	}
	fp=fopen(outfile2,"w");
	fprintf(fp,"#t rho abs(m) m\n");
	fclose(fp);

	sprintf(readfile1,"acceptprobK%dL%dM%5.4lf",K,L,mu);
	fp=fopen(readfile1,"r");
	if(fp==NULL)
	{
		printf("The FILE [%s] DOES NOT EXIST\n",readfile1);
		exit(0);
	}
	while(fscanf(fp,"%d%lf",&i,&x)!=EOF)
		acceptance[i]=x;
	fclose(fp);
	if(i!=L)
	{
		printf("Error in FILE [%s]\n",readfile1);
		exit(0);
	}

	sprintf(readfile2,"periodicprobK%dL%dM%5.4lf",K,L,mu);
	fp=fopen(readfile2,"r");
	if(fp==NULL)
	{
		printf("The FILE [%s] DOES NOT EXIST\n",readfile2);
		exit(0);
	}
	while(fscanf(fp,"%d%lf",&i,&x)!=EOF)
		periodic[i]=x;
	fclose(fp);
	if(i!=K-1)
	{
		printf("Error in FILE [%s]\n",readfile2);
		exit(0);
	}

}

void remove_hor()
{
	/* removes all horizontal kmers */
	int i;

	for(i=0;i<lasth;i++)
		lat[horizontal[i]]=0;
	lasth=0;
}

void remove_ver()
{
	/* removes all vertical kmers */
	int i;

	for(i=0;i<lastv;i++)
		lat[vertical[i]]=0;
	lastv=0;
}

void deposit_hor(int i)
{
	/*puts a horizontal kmer with head at i*/

	int j;
	
	for(j=1;j<=K;j++)
	{
		lat[i]=1;
		horizontal[lasth]=i;
		lasth++;
		i=rn[i];
	}
}

void deposit_ver(int i)
{
	/*puts a vertical kmer with head at i*/

	int j;
	
	for(j=1;j<=K;j++)
	{
		lat[i]=2;
		vertical[lastv]=i;
		lastv++;
		i=tn[i];
	}
}

void fill_periodic_hor(int j)
{
	/*Given an empty horizontal line, the function checks the
	 * occupation probability of first K-1 sites */

	int i,k;

	starth=j;finalh=ln[j];
	for(i=0;i<K-1;i++)
	{
		if(genrand_real3() < periodic[i])
		{
			for(k=1;k<=K;k++)
			{
				lat[starth]=1;
				horizontal[lasth]=starth;
				lasth++;
				starth=rn[starth];
			}
			return;
		}
		else
		{
			finalh=starth;
			starth=rn[starth];
		}
	}
	return;
}

void fill_periodic_ver(int j)
{
	/*Given an empty vertical line, the function checks the
	 * occupation probability of first K-1 sites */

	int i,k;

	startv=j;finalv=bn[j];
	for(i=0;i<K-1;i++)
	{
		if(genrand_real3() < periodic[i])
		{
			for(k=1;k<=K;k++)
			{
				lat[startv]=2;
				vertical[lastv]=startv;
				lastv++;
				startv=tn[startv];
			}
			return;
		}
		else
		{
			finalv=startv;
			startv=tn[startv];
		}
	}
	return;
}

int find_starthfinalh(int row)
{
	/*in row 'row' finds out starth and finalh*/

	int i;

	starth=row*L; finalh=row*L;
	if(lat[starth]!=0)
	{
		while(lat[starth]!=0)
		{
			starth=rn[starth];
			if(starth==row*L)
				return 0;
		}
		while(lat[finalh]!=0)
			finalh=ln[finalh];
		return 1;
	}

	while(lat[ln[starth]]==0)
	{
		starth=ln[starth];
		if(starth==row*L)
		{
			fill_periodic_hor(starth);
			return 1;
		}
	}
	finalh=ln[starth];
	while(lat[finalh]!=0)
		finalh=ln[finalh];
	return 1;
}

void fill_row(int row)
{
	/*fills row 'row' with horizontal kmers*/

	int i,len;

	if(find_starthfinalh(row)==0)
		return;
	do
	{
		endh=starth;len=1;
		while(lat[rn[endh]]==0)
		{
			endh=rn[endh];
			len++;
			if(endh==finalh)
				break;
		}
		while(len>=K)
		{
			if(genrand_real3() < acceptance[len])
			{
				deposit_hor(starth);
				for(i=0;i<K;i++)
					starth=rn[starth];
				len=len-K;
			}
			else
			{
				starth=rn[starth];
				len--;
			}
		}
		if(endh==finalh)
			return;
		starth=rn[endh];
		while(lat[starth]!=0)
			starth=rn[starth];
	}
	while(endh!=finalh);
}

int find_startvfinalv(int col)
{
	/*in row 'row' finds out startv and finalv*/

	int i;

	startv=col; finalv=col;
	if(lat[startv]!=0)
	{
		while(lat[startv]!=0)
		{
			startv=tn[startv];
			if(startv==col)
				return 0;
		}
		while(lat[finalv]!=0)
			finalv=bn[finalv];
		return 1;
	}

	while(lat[bn[startv]]==0)
	{
		startv=bn[startv];
		if(startv==col)
		{
			fill_periodic_ver(startv);
			return 1;
		}
	}
	finalv=bn[startv];
	while(lat[finalv]!=0)
		finalv=bn[finalv];
	return 1;
}

void fill_col(int col)
{
	/*fills col 'col' with horizontal kmers*/

	int i,len;

	if(find_startvfinalv(col)==0)
		return;
	do
	{
		endv=startv;len=1;
		while(lat[tn[endv]]==0)
		{
			endv=tn[endv];
			len++;
			if(endv==finalv)
				break;
		}
		while(len>=K)
		{
			if(genrand_real3() < acceptance[len])
			{
				deposit_ver(startv);
				for(i=0;i<K;i++)
					startv=tn[startv];
				len=len-K;
			}
			else
			{
				startv=tn[startv];
				len--;
			}
		}
		if(endv==finalv)
			return;
		startv=tn[endv];
		while(lat[startv]!=0)
			startv=tn[startv];
	}
	while(endv!=finalv);
}

void evolve()
{
	/*FILE *fp;
	fp=fopen("config_L42.dat","w");*/
	int row,col,i,id,u,v,dest,source;
	/*if(taskid == MASTER)
	{
		for(dest=1;dest<numtasks;dest++)
			MPI_Send(&lat[L*BLOCK_LOW(dest,numtasks,L)], L*BLOCK_SIZE(dest,numtasks,L), MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
	}
	else
	MPI_Recv(&lat[L*BLOCK_LOW(taskid,numtasks,L)], L*BLOCK_SIZE(taskid,numtasks,L), MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);*/

	MPI_Scatterv(lat,scounts,displs,MPI_INT,&lat[L*BLOCK_LOW(taskid,numtasks,L)],L*BLOCK_SIZE(taskid,numtasks,L),MPI_INT,MASTER,MPI_COMM_WORLD);
	remove_hor();
	for(row=BLOCK_LOW(taskid,numtasks,L);row<=BLOCK_HIGH(taskid,numtasks,L);row++)
		fill_row(row);
	MPI_Gatherv(&lat[L*BLOCK_LOW(taskid,numtasks,L)],L*BLOCK_SIZE(taskid,numtasks,L),MPI_INT,lat,scounts,displs,MPI_INT,MASTER,MPI_COMM_WORLD);
	/*if(taskid != MASTER)
		MPI_Send(&lat[L*BLOCK_LOW(taskid,numtasks,L)], L*BLOCK_SIZE(taskid,numtasks,L), MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD);
	else 
	{		
		for(source=1;source<numtasks;source++)		
		MPI_Recv(&lat[L*BLOCK_LOW(source,numtasks,L)], L*BLOCK_SIZE(source,numtasks,L), MPI_INT, source, FROM_WORKER, MPI_COMM_WORLD, &status);
	}*/	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(taskid == MASTER)
	{
		for(dest=1;dest<numtasks;dest++)
			MPI_Send(&lat[BLOCK_LOW(dest,numtasks,L)], 1, columntype, dest, FROM_MASTER, MPI_COMM_WORLD);
	}
	else
		MPI_Recv(&lat[BLOCK_LOW(taskid,numtasks,L)], 1, columntype, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
	remove_ver();
	for(col=BLOCK_LOW(taskid,numtasks,L);col<=BLOCK_HIGH(taskid,numtasks,L);col++)
		fill_col(col);
	if(taskid != MASTER)
		MPI_Send(&lat[BLOCK_LOW(taskid,numtasks,L)], 1, columntype, MASTER, FROM_WORKER, MPI_COMM_WORLD);
	else 
	{		
		for(source=1;source<numtasks;source++)		
			MPI_Recv(&lat[BLOCK_LOW(source,numtasks,L)], 1, columntype, source, FROM_WORKER, MPI_COMM_WORLD, &status);
	}

	
	/*if(taskid == 1)
	{	
		for(u=L-1;u>=0;u--)
		{
			for(v=0;v<L;v++)	
				fprintf(fp,"%d ",lat[u*L+v]);
			fprintf(fp,"\n");
		}
	}
	fclose(fp);*/
}


void take_reading(int j)
{
	int i;
	double tmp,sum;

	tmp=1.0*(lasth0-lastv0)/(1.0*N);
	meanrho[j]=meanrho[j]+1.0*(lasth0+lastv0)/(1.0*N);
	rho2[j]=rho2[j]+1.0*(lasth0+lastv0)*(lasth0+lastv0)/(1.0*N*N);
	meanabs[j]=meanabs[j]+fabs(tmp);
	meanm1[j]=meanm1[j]+tmp;
	meanm2[j]=meanm2[j]+tmp*tmp;
	meanm4[j]=meanm4[j]+tmp*tmp*tmp*tmp;
	i=floor((tmp+1.0)/BINSIZE);
	if(lasth0-lastv0==-N)
		i=0;
	if(lasth0-lastv0==N)
		i=BINS-1;
	if((i<BINS) && (i>=0))
		count[i]++;
}

void output2(int t)
{
	double tmp;
	FILE *fp;

	tmp=1.0*(lasth0-lastv0)/(1.0*N);
	fp=fopen(outfile2,"a");
	fprintf(fp,"%d %e %e %e\n",t,1.0*(lastv0+lasth0)/(1.0*N),fabs(tmp),tmp);
	fclose(fp);
}

void output1(int ave,int window)
{
	int i;
	double sum1,sum2,tmp,tmp1,tmp2,avgrho;
	FILE *fp;

	tmp=1.0/(1.0*ave);
	tmp1=1.0/(window*1.0);
	tmp2=sqrt(tmp1);
	for(i=0;i<window;i++)
	{
		Rho[i]=meanrho[i]*tmp;
		Rho_2[i]=rho2[i]*tmp;
	}
	for(i=0;i<window;i++)
	{
		modQ[i]=meanabs[i]*tmp/Rho[i];
		Q_1[i]=meanm1[i]*tmp/Rho[i]; 
		Q_2[i]=meanm2[i]*tmp/pow(Rho[i],2.0); 
		Q_4[i]=meanm4[i]*tmp/pow(Rho[i],4.0);
	}
	for(i=0;i<window;i++)
	{
		Rho_2[i]=(Rho_2[i]-Rho[i]*Rho[i])*L*L; 
		Q_4[i]=1.0-Q_4[i]/Q_2[i]/Q_2[i]/3.0;
		Q_2[i]=(Q_2[i]-modQ[i]*modQ[i])*L*L; 
	}
	fp=fopen(outfile1,"w");
	fprintf(fp,"%e	%e",1.0*window*ave,mu);

	sum1=0;sum2=0;
	for(i=0;i<window;i++)
	{
		sum1=sum1+Rho[i];
		sum2=sum2+Rho[i]*Rho[i];
	}
	sum1=sum1*tmp1;sum2=sum2*tmp1;
	fprintf(fp," %e %e",sum1,sqrt(sum2-sum1*sum1)*tmp2);
	sum1=0;sum2=0;
	for(i=0;i<window;i++)
	{
		sum1=sum1+Rho_2[i];
		sum2=sum2+Rho_2[i]*Rho_2[i];
	}
	sum1=sum1*tmp1;sum2=sum2*tmp1;
	fprintf(fp," %e %e",sum1,sqrt(sum2-sum1*sum1)*tmp2);

	sum1=0;sum2=0;
	for(i=0;i<window;i++)
	{
		sum1=sum1+modQ[i];
		sum2=sum2+modQ[i]*modQ[i];
	}
	sum1=sum1*tmp1;sum2=sum2*tmp1;
	fprintf(fp," %e %e",sum1,sqrt(sum2-sum1*sum1)*tmp2);

	sum1=0;sum2=0;
	for(i=0;i<window;i++)
	{
		sum1=sum1+Q_1[i];
		sum2=sum2+Q_1[i]*Q_1[i];
	}
	sum1=sum1*tmp1;sum2=sum2*tmp1;
	fprintf(fp," %e %e",sum1,sqrt(sum2-sum1*sum1)*tmp2);

	sum1=0;sum2=0;
	for(i=0;i<window;i++)
	{
		sum1=sum1+Q_2[i];
		sum2=sum2+Q_2[i]*Q_2[i];
	}
	sum1=sum1*tmp1;sum2=sum2*tmp1;
	fprintf(fp," %e %e",sum1,sqrt(sum2-sum1*sum1)*tmp2);

	sum1=0;sum2=0;
	for(i=0;i<window;i++)
	{
		sum1=sum1+Q_4[i];
		sum2=sum2+Q_4[i]*Q_4[i];
	}
	sum1=sum1*tmp1;sum2=sum2*tmp1;
	fprintf(fp," %e %e\n",sum1,sqrt(sum2-sum1*sum1)*tmp2);

	fclose(fp);


	fp=fopen(outfile3,"w");
	fprintf(fp,"#m P(m)\n");
	for(i=0;i<BINS;i++)
	{
		if(count[i]!=0)
			fprintf(fp,"%e %e\n",mass[i],count[i]/factor[i]/ave/window);
	}
	fclose(fp);

}

void chk_config()
{
	int i,j;
	j=0;
	for (i=0;i<N;i++)
	{
		if(lat[i]!=0)
			j++;
		if((lat[i]>= DIR)||(lat[i]<0))
		{
			printf("ERROR IN STATE\n");
			sleep(3);
		}
	}
	if(j!=(lastv0+lasth0))
	{
		printf("ERROR IN BOOK KEEPING\n");
		sleep(3);
	}
}

int main (int argc, char *argv[])
{
	int source,dest,mtype,ave,extra,offset1,offset2,i,j,g,rc,m,n,u,k;           

	MPI_Init(&argc,&argv);
	MPI_Barrier(MPI_COMM_WORLD);
	/*elapsed_time=-MPI_Wtime();*/
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

	lat = (int *) malloc (N * sizeof(int));
	ln = (int *) malloc (N * sizeof(int));
	rn = (int *) malloc (N * sizeof(int));
	tn = (int *) malloc (N * sizeof(int));
	bn = (int *) malloc (N * sizeof(int));
	horizontal = (int *) malloc (N * sizeof(int));
	vertical = (int *) malloc (N * sizeof(int));
	displs = (int *) malloc(numtasks * sizeof(int));
	scounts = (int *) malloc(numtasks * sizeof(int));
	
	take_input();
	neighbour();
	/*lat_initialization();*/
	if (taskid == 0)
	{		
		lat_initialization();
		mu=5.62;
		initialize();
		lat_init();
	}
	MPI_Bcast (&mu, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Bcast (&acceptance[0], L+1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Bcast (&periodic[0], K, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	
	MPI_Type_vector(L, BLOCK_SIZE(taskid,numtasks,L), L, MPI_INT, &columntype);
	MPI_Type_commit(&columntype);

	for(u=0;u<numtasks;u++)
	{
		displs[u] = (u*L*BLOCK_SIZE(u,numtasks,L));
		scounts[u]= L*BLOCK_SIZE(u,numtasks,L);
	}

	for(i=0;i<T;i++)
	{
		evolve();
		MPI_Reduce(&lasth, &lasth0, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
		MPI_Reduce(&lastv, &lastv0, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
		if(taskid == MASTER)
		{	if(i % GAP2 ==0)							
				output2(i);
		}
	}

	for(j=0;j<BLOCKS;j++)
	{
		for(i=0;i<AVE;i++)
		{
			evolve();
			MPI_Reduce(&lasth, &lasth0, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
			MPI_Reduce(&lastv, &lastv0, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
			if(taskid == MASTER)
			{
				if(i % GAP2 ==0)
					output2(T+j*AVE+i);
				take_reading(j);
			}
		}
		if(taskid == MASTER)
		{
			g=j+1;
			output1(AVE,g);
		}
	}

	MPI_Finalize();
	
}


