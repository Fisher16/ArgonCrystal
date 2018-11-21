#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include "argon.hpp"

//TODO: rooll up some leftover x y z code duplicates



std::vector<vec> init_pos(int N,std::string filename){
	std::vector<vec> r(N);

	std::ofstream init_cond;

	init_cond.open(filename);
	init_cond<<N<<"\n\n";

	for(int i0=0;i0<n;i0++)
		for(int i1=0;i1<n;i1++)
			for(int i2=0;i2<n;i2++){
				int i=i0+i1*n+i2*n*n;
				//From (5) pdf
				r[i]=(i0-(n-1)/2)*b0+(i1-(n-1)/2)*b1+(i2-(n-1)/2)*b2;
				init_cond<<"Ar "<<r[i].x<<" "<<r[i].y<<" "<<r[i].z<<"\n";
			}
	init_cond.close();
	return r;
}
void append(std::string filename,std::vector<vec> r_list){
	std::ofstream file;
	file.open(filename,std::ios_base::app);
	file<<N<<"\n\n";
		for(auto r : r_list)file<<"Ar "<<r.x<<" "<<r.y<<" "<<r.z<<"\n";
	file.close();
}



std::vector<vec> init_momentum(int N,double T0,std::string filename){
	std::vector<vec> p(N);

	std::ofstream init_cond;
	init_cond.open(filename);

	std::vector<vec> Ek_list(N);
	//TODO: check if short int is faster
	std::vector<sign> sign_list(N);

	//lambda in lmabda
	auto gen_Ek =[T0](){
		double lambda=(double)rand() / (RAND_MAX);
		return -0.5*k*T0*log(lambda);
	};

	//get distribution
	vec E_theoretical{0.5*k*T0,0.5*k*T0,0.5*k*T0};
	vec E_dist{0,0,0};
	for(int i=0;i<N;i++){
		Ek_list[i].x=gen_Ek();
		sign_list[i].x=2*(rand()%2)-1;

		Ek_list[i].y=gen_Ek();
		sign_list[i].y=2*(rand()%2)-1;

		Ek_list[i].z=gen_Ek();
		sign_list[i].z=2*(rand()%2)-1;
		E_dist=E_dist+Ek_list[i];
	}
	std::cout<<"Distribution\n";
	print(E_theoretical);
	E_dist=((double)1/N)*E_dist;
	print(E_dist);
	vec E_diff=E_dist-E_theoretical;
	print(E_diff);

	E_dist={0,0,0};
	for(int i=0;i<N;i++){
		Ek_list[i]=Ek_list[i]-E_diff;
		Ek_list[i].x=Ek_list[i].x<0?0.0:Ek_list[i].x;
		Ek_list[i].y=Ek_list[i].y<0?0.0:Ek_list[i].y;
		Ek_list[i].z=Ek_list[i].z<0?0.0:Ek_list[i].z;

		E_dist=E_dist+Ek_list[i];
	}

	E_dist=((double)1/N)*E_dist;
	print(E_dist);
	//calc momentum
	vec P{0,0,0};
	for(int j=0;j<N;j++){
		p[j].x=sign_list[j].x*sqrt(2*m*Ek_list[j].x);
		p[j].y=sign_list[j].y*sqrt(2*m*Ek_list[j].y);
		p[j].z=sign_list[j].z*sqrt(2*m*Ek_list[j].z);
		P=P+p[j];
	}

	//normalize momentum
	//std::cout<<P.x<<" "<<P.y<<" "<<P.z<<"\n";
	P=((double)1/N)*P;
	for(int i=0;i<N;i++){
		p[i]=p[i]-P;
		init_cond<<p[i].x<<" "<<p[i].y<<" "<<p[i].z<<"\n";
	}

	init_cond.close();
	return p;
}

vec VDW(vec ri, vec rj){
	double rij=dist(ri,rj);
	return (12*epsilon/rij/rij*(pow(R/rij,12)-pow(R/rij,6)))*(ri-rj);
}


vec BALL(vec ri){
	double r_len=len(ri);
	if(r_len>=L)
		return (f*(L-r_len)/r_len)*ri;
	else
		return {0.00,0.00,0.00};

}

double Ek(vec p){
	return (0.5/m)*sq(p);
}

double Ek_tot(std::vector<vec> p_list){
	double Ek_tot=0;
	for(auto p : p_list)Ek_tot+=Ek(p);
	return Ek_tot;
}

double calc_T(double Ek_tot){
	return 2*Ek_tot/(3*N*k);

}

double Ep_wdv(double rij){
	return epsilon*(pow(R/rij,12)-2*pow(R/rij,6));
}

double Ep_ball(double r_len){
	return r_len<L?0:0.5*f*(r_len-L)*(r_len-L);
}

double Ep_tot(std::vector<vec> r_list){
	double Ep=0;
	for(int i=0;i<N-1;i++){
		for(int j=i;j<N;j++){
			if(i!=j){
				Ep+=Ep_wdv(dist(r_list[i],r_list[j]));
			}
		}
	}
	//std::cout<<"E vdw"<<Ep<<"\n";
	for(int i=0;i<N;i++)Ep+=Ep_ball(len(r_list[i]));

	return Ep;
}



double E_tot(std::vector<vec> r_list,std::vector<vec> p_list){
	//std::cout<<Ek_tot(p_list)<<"\n"<<Ep_tot(r_list)<<"\n";
	return Ep_tot(r_list)+Ek_tot(p_list);
}

std::vector<vec> calc_force(std::vector<vec> r_list){
	int const N=r_list.size();
	std::vector<vec> F(N);
	for(int j=0;j<N;j++)F[j]={0,0,0};

	for(int i=0;i<N-1;i++){
		for(int j=i;j<N;j++){
			if(i!=j){
				vec force=VDW(r_list[i],r_list[j]);
				F[i]=(F[i]+force);
				F[j]=(F[j]-force);
			}
		}
		F[i]=F[i]+BALL(r_list[i]);
	}

	F[N-1]=F[N-1]+BALL(r_list[N-1]);
	return F;
}


std::vector<vec> time_step(double const dt,std::vector<vec> &r_list, std::vector<vec> &p_list, std::vector<vec> &F_list){
	//TODO: check if faster multiply each elem and then add
	for(int i=0;i<N;i++){
		p_list[i]=p_list[i]+((dt/2)*F_list[i]);
		r_list[i]=r_list[i]+((dt/m)*p_list[i]);
	}
	std::vector<vec> Fdt=calc_force(r_list);
	for(int i=0;i<N;i++)
		p_list[i]=p_list[i]+((dt/2)*Fdt[i]);

	return Fdt;
}

int main(){
	srand(12345);
	double T0=1000;
	//double t=0.0;
	double const dt=0.0001;
	while(T0<1100){
		auto T0str = std::to_string(T0);
		auto Nstr = std::to_string(N);
		std::vector<vec> r_list=init_pos(N,"data/init_pos.txt");
		std::vector<vec> p_list=init_momentum(N,T0,"data/init_momentum.txt");
		//if(abs(T0)<0.0001)for(int i=0;i<N;i++)p_list[i]={0,0,0};
		std::vector<vec> F_list=calc_force(r_list);

		std::cout<<T0<<"\n";
		std::cout<<E_tot(r_list,p_list)<<"\n";

		double sim_time=1;

		std::vector<vec> F_1(N);
		std::vector<vec> F_2(F_list.begin(),F_list.end());

		int i=0;
		std::string simout="data/sim"+Nstr+"_"+T0str+".txt";
		std::string ETout="data/E_T_"+Nstr+"_"+T0str+".txt";
		//clearout previous
		std::ofstream f(simout);
		f.close();
		std::ofstream ff(ETout);
		ff.close();
		std::ofstream fET(ETout,std::ios_base::app);
		//TODO: made simulate function
		for(double t=0;t<sim_time;t+=2*dt){
			F_1=time_step(dt,r_list,p_list,F_2);
			F_2=time_step(dt,r_list,p_list,F_1);
			if(i%100==0){
				append(simout,r_list);
				i=0;
			}
			if(i%10==0){
				fET<<t<<" "<<E_tot(r_list,p_list)<<" "<<Ek_tot(p_list)<<" "<<Ep_tot(r_list)<<" ";
				fET<<calc_T(Ek_tot(p_list))<<"\n";
			}
			i+=2;

		}

		std::cout<<E_tot(r_list,p_list)<<"\nTemp ";
		std::cout<<calc_T(Ek_tot(p_list))<<"\n";
		fET.close();
		f.close();
		T0+=100;
	}
	return 0;
}
