//DATA STRUCTURES
struct vec{
	double x;
	double y;
	double z;
};

struct sign{
  int8_t x;
  int8_t y;
  int8_t z;
};

vec operator+(vec a, vec b){
	return vec{a.x+b.x,a.y+b.y,a.z+b.z};
}
vec operator-(vec a, vec b){
	return vec{a.x-b.x,a.y-b.y,a.z-b.z};
}

vec operator*(double val, vec a){
	return vec{val*a.x,val*a.y,val*a.z};
}
double sqr(vec a){
	return a.x*a.x+a.y*a.y+a.z*a.z;
}
double dist(vec a, vec b){
	return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

double len(vec a){
	return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

double sq(vec a){
	return a.x*a.x+a.y*a.y+a.z*a.z;
}

void print(vec a){
	std::cout<<"X: "<<a.x<<" Y: "<<a.y<<" Z: "<<a.z<<"\n";
}



//PHYSICS CONSTANTS
double const a=0.38;
double const k=8.31/1000;
double R=a;
double epsilon=1;

double const m=1;


//SIMULATION PARAMETERS
int const n=5;
int const N=n*n*n;
//Preseure sim sphere
double const L=1.513158*a*(n-1);
double const f=10000.0;

//edges of elementar cell
//TODO: use more precies values if needed

vec b0{a,0,0};
vec b1{0.5*a,0.86602540378*a,0};
vec b2{0.5*a,0.28867513459*a,0.81649658092*a};

/*
vec b0{a,0,0};
vec b1{0,a,0};
vec b2{0,0,a};
*/
