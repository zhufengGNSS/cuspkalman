//fieldmap.cxx

#include "fieldmap.h"

double Er_min,Er_max,Er_step,
       Ez_min,Ez_max,Ez_step;

double Br_min,Br_max,Br_step,
       Bz_min,Bz_max,Bz_step;       

/* 
 Adjust the settings of these matrices: rows,columns
 For example for the WITCH experiment: E:(1000,1000), B:(701,26)
 For the ASACUSA experiment: E:(1000,1000), B:(500,500)
*/

const double Enrows = 1000;
const double Encols = 1000;
const double Bnrows = 1000;
const double Bncols = 1000;
Matrix<pair<double, double> > * Ematrix = new Matrix<std::pair<double, double> >(Enrows,Encols); //(Ez,Er)
Matrix<pair<double, double> > * Bmatrix = new Matrix<std::pair<double, double> >(Bnrows,Bncols);    //(Bz,Br)


//pool things
double r,dzr;
double z_index,r_index;
double z1,z2,r1,r2;
int z1_index,z2_index,r1_index,r2_index;

pair<double, double> B12,B22,B11,B21, Bzr;
pair<double, double> E12,E22,E11,E21, Ezr;

vector<double> electric_field(3);
vector<double> magnetic_field(3);
    
    
/*
int round(&double num) {
  char sign = static_cast<char>(num/fabs(num));
  return static_cast<int>(sign* (fabs(num)+ 0.5));
}    */

vector<double> values;
 


double get_Ez_min(){return Ez_min;}
double get_Ez_max(){return Ez_max;}
double get_Ez_step(){return Ez_step;}
double get_Er_min(){return Er_min;}
double get_Er_max(){return Er_max;}
double ger_Er_step(){return Er_step;}

double get_Bz_min(){return Bz_min;}
double get_Bz_max(){return Bz_max;}
double get_Bz_step(){return Bz_step;}
double get_Br_min(){return Br_min;}
double get_Br_max(){return Br_max;}
double ger_Br_step(){return Br_step;}
 

bool IsElement(double _value){

  bool returnvalue = false;
  
  vector<int>::reverse_iterator rit;
  for ( unsigned int index= 0 ; index < values.size(); index++ ){
    if ( values[index] == _value) returnvalue = true;
  }
  return returnvalue;
}  

Matrix<pair<double, double> > * ReadElectricFieldAscii(const char * _Erz_filename){
     //Reads from a file constructed as r z Er Ez separated by tabs or spaces.
     //For each r stepped, all z-stepped are being walked through.
     //steps over z-values first
     //should start with z-values from low to high values    
    ifstream Erz_str;
    
    Erz_str.open(_Erz_filename); 
    if (!Erz_str){cout<<"Error in opening "<<_Erz_filename<<". The program will be stopped."<<endl; exit(1);}
     
    double r,z,Ez,Er;//,Emod;
    //     double r2,z2;
     int r_stepper = 0, z_stepper = 0;
     char firstline[256]; 
     char comment= '%';
     while (Erz_str.peek() == comment){
     	Erz_str.getline(firstline, 256);
     }
     
     //read file minus first whiteline.
     bool valuesset = false;

     while ( !Erz_str.eof() && z_stepper != Enrows && r_stepper != Encols ) { // keep reading until end-of-file 
        Erz_str>>r>>z>>Er>>Ez;
        //Br_min,max en voor z ook bepalen he.
        if(!valuesset){valuesset=true;Er_min=Er_max=r;Ez_min=Ez_max=z;}
        if(r<Er_min){Er_min=r;}if(r>Er_max){Er_max=r;}if(z<Ez_min){Ez_min=z;} if(z>Ez_max){Ez_max=z;}        
        (*Ematrix)(z_stepper,r_stepper).first=Ez;
	(*Ematrix)(z_stepper,r_stepper).second=Er;
        z_stepper++;
        if (z_stepper == Enrows){z_stepper =0;r_stepper++;}                     
     }
     
     
     
     //set the r and z_step
     Er_step=(Er_max-Er_min)/(Encols-1);
     Ez_step=(Ez_max-Ez_min)/(Enrows-1);
     
     cout<<"\t R: "<<Er_min<<" "<<Er_max<<" "<<Er_step<<endl;
     cout<<"\t z: "<<Ez_min<<" "<<Ez_max<<" "<<Ez_step<<endl;
     /*
     for (unsigned i=0; i < Ematrix.nrows(); i++){
         for (unsigned j=0; j < Ematrix.ncols();j++){
             cout<<Ez_min+i*Ez_step<<" \t"<<Er_min+j*Er_step<<" \t"<<Ematrix(i,j).first<<endl;
             }
             cin.get();
     }    
    */
    cout<<"file: "<<_Erz_filename<<" read\n"; 
    Erz_str.close();

    return Ematrix;
}    
    
    
/*read function*/    
Matrix<pair<double, double> > * ReadElectricFieldAscii(const char *_Er_filename,const char *_Ez_filename){
     //Supposed it that the 2 files have the same number of indexes!!!
    ifstream Er_str, Ez_str;
    
    Er_str.open(_Er_filename); 
    Ez_str.open(_Ez_filename);
    if (!Er_str){cout<<"Error in opening "<<_Er_filename<<". The program will be stopped."<<endl; exit(1);}
    if (!Er_str){cout<<"Error in opening "<<_Er_filename<<". The program will be stopped."<<endl; exit(1);}
     
    double r,z,Ez,Er;//,Emod;
    //     double r2,z2;
     int r_stepper = 0, z_stepper = 0;
     char firstline[256]; 
     char comment= '%';
     while (Er_str.peek() == comment){
     	Er_str.getline(firstline, 256);
     }
     while (Ez_str.peek() == comment){
     	Ez_str.getline(firstline, 256);
     }
     
     //read file minus first whiteline.
     bool valuesset = false;


     while ( !Er_str.eof() && z_stepper != Enrows && r_stepper != Encols ) { // keep reading until end-of-file 
        Er_str>>r>>z>>Er;
        Ez_str>>r>>z>>Ez;
        //Br_min,max en voor z ook bepalen he.
        if(!valuesset){valuesset=true;Er_min=Er_max=r;Ez_min=Ez_max=z;}
        if(r<Er_min){Er_min=r;}if(r>Er_max){Er_max=r;}if(z<Ez_min){Ez_min=z;} if(z>Ez_max){Ez_max=z;}        
        (*Ematrix)(z_stepper,r_stepper).first=Ez;
        (*Ematrix)(z_stepper,r_stepper).second=Er;
        z_stepper++;
        if (z_stepper == Enrows){z_stepper =0;r_stepper++;}                     
     }
     
     //set the r and z_step
     Er_step=(Er_max-Er_min)/(Encols-1);
     Ez_step=(Ez_max-Ez_min)/(Enrows-1);
     
     cout<<"\t R: "<<Er_min<<" "<<Er_max<<" "<<Er_step<<endl;
     cout<<"\t z: "<<Ez_min<<" "<<Ez_max<<" "<<Ez_step<<endl;
     /*
     for (unsigned i=0; i < Ematrix.nrows(); i++){
         for (unsigned j=0; j < Ematrix.ncols();j++){
             cout<<Ez_min+i*Ez_step<<" \t"<<Er_min+j*Er_step<<" \t"<<Ematrix(i,j).first<<endl;
             }
             cin.get();
     }    
    */
    cout<<"file: "<<_Er_filename<<" read\n";
    cout<<"file: "<<_Ez_filename<<" read\n";    
    Er_str.close(); 
    Ez_str.close();

    return Ematrix;

}

Matrix<pair<double, double> > * ReadMagneticFieldAscii(const char *_B_filename){
  //Reads from a file constructed as r z Br Bz separated by tabs or spaces.
  //For each r stepped, all z-stepped are being walked through.
  //steps over z-values first
  //should start with z-values from low to high values
  ifstream B_str;
  B_str.open(_B_filename); 
  if (!B_str){cout<<"Error in opening "<<_B_filename<<". The program will be stopped."<<endl; exit(1);}
  double r,z,Bz,Br;
  int r_stepper = 0, z_stepper = 0;
  
  
  bool valuesset = false;
  
  char firstline[256]; 
  char comment= '%';
  while (B_str.peek() == comment){
    B_str.getline(firstline, 256);
    //	cout<<firstline<<endl;
  }
  
  int i =0;
  while ( !B_str.eof() && z_stepper != Bnrows && r_stepper != Bncols ) { // keep reading until end-of-file 
    B_str>>r>>z>>Br>>Bz;
    //cout<<i<<" "<<z<<" "<<r<<" "<<Bz<<" "<<Br<<" "<<endl;cin.get();
    
    //Br_min,max en voor z ook bepalen he.
    if(!valuesset){valuesset=true;Br_min=Br_max=r;Bz_min=Bz_max=z;}
    if(r<Br_min){Br_min=r;}if(r>Br_max){Br_max=r;}if(z<Bz_min){Bz_min=z;} if(z>Bz_max){Bz_max=z;}        
    (*Bmatrix)(z_stepper,r_stepper).first=Bz;
    (*Bmatrix)(z_stepper,r_stepper).second=Br;
    z_stepper++;
    if (z_stepper == Bnrows){;z_stepper =0;r_stepper++;}   
    i++;                  
  }
  
  //set the r and z_step
  Br_step=(Br_max-Br_min)/(Bncols-1);
  Bz_step=(Bz_max-Bz_min)/(Bnrows-1);
  
  cout<<"\t R: "<<Br_min<<" "<<Br_max<<" "<<Br_step<<endl;
  cout<<"\t z: "<<Bz_min<<" "<<Bz_max<<" "<<Bz_step<<endl;
  /*
    for (unsigned i=0; i < Bmatrix.nrows(); i++){
    for (unsigned j=0; j < Bmatrix.ncols();j++){
    cout<<Bz_min+i*Bz_step<<" \t"<<Br_min+j*Br_step<<" \t"<<Bmatrix(i,j).first<<endl;
    }
    cin.get();
    }    
  */
  cout<<"file: "<<_B_filename<<" read\n";
  B_str.close(); 
  return Bmatrix;
}


  
/*Interpolate functions*/  
vector<double> getElectrAscii(const double& _x,const double& _y,const double& _z){

     //Interpolate potentials on r,z direction and calculate stuff
     r=sqrt(_x*_x+_y*_y);if(r == 0)r=1e-10; //otherwise division by zero!     
     //check of ze wel in of bound zitten he!
     if(_x>Er_max||_y>Er_max||r>Er_max ||_z>Ez_max||_z<(Ez_min+Er_step)){
           cout<<"E particle out of bound: "<<_x<<" "<<_y<<" "<<_z<<endl;}
     //scale both particles
      z_index= _z-Ez_min; 
      r_index= r-Er_min;


      r1_index=fabs(floor(r_index/Er_step));r2_index=r1_index+1;//ceil(r_index/Er_step); //the next one he: (r+0.5*r_step)/r_step
      z1_index=floor(z_index/Ez_step);z2_index=z1_index+1;//ceil(z_index/Ez_step);

      //matrix:first row(z), then column(r)
      
      E12=(*Ematrix)(z2_index,r1_index);E22=(*Ematrix)(z2_index,r2_index);
      E11=(*Ematrix)(z1_index,r1_index);E21=(*Ematrix)(z1_index,r2_index);
      
      r1=r1_index*Er_step+Er_min;r2=r2_index*Er_step+Er_min;
      z1=z1_index*Ez_step+Ez_min;z2=z2_index*Ez_step+Ez_min;
      
   
      //interpolate Potential!  http://en.wikipedia.org/wiki/Bilinear_interpolation
      dzr=1.0/((r2-r1)*(z2-z1));
      Ezr.first = E11.first*(r2-r)*(z2-_z)*dzr
                 +E21.first*(r-r1)*(z2-_z)*dzr
                 +E12.first*(r2-r)*(_z-z1)*dzr
                 +E22.first*(r-r1)*(_z-z1)*dzr;
      Ezr.second = E11.second*(r2-r)*(z2-_z)*dzr
                 +E21.second*(r-r1)*(z2-_z)*dzr
                 +E12.second*(r2-r)*(_z-z1)*dzr
                 +E22.second*(r-r1)*(_z-z1)*dzr;
      //if nan, solve it this way he!           
      if(Ezr.second != Ezr.second){ 
            Ezr.second = 0.0;
            cout<<"Efoutje: "<<dzr<<endl;
            cout<<"x: "<<_x<<endl;
            cout<<"y: "<<_y<<endl;
            cout<<"r: "<<r<<endl;
            cout<<"z: "<<_z<<endl; 
      cout<<r1_index<<"\t"<<r2_index<<endl;
      cout<<z1_index<<"\t"<<z2_index<<endl;            
            cin.get();
            electric_field[0]=0.0;//E11.second;
            electric_field[1]=0.0;
      }
      else{
           /*if (fabs(_x) < 0.0001)*/ electric_field[0]=Ezr.second*(_x/r); //else electric_field[0]=0.0;
           /*if (fabs(_y) < 0.0001)*/ electric_field[1]=Ezr.second*(_y/r); //else electric_field[0]=0.0;
      }
      
      if(Ezr.first != Ezr.first){ 
           Ezr.first=0.0;
           electric_field[2]=0.0;                                  
      }
      else{
           electric_field[2]=Ezr.first;           
     }
        //cout<<"asked for E with z="<<_z<<" r="<<r<<"  ->Ez="<<Ezr.first<<" // Er="<<Ezr.second<<" V/m"<<endl;
     return electric_field;                
}     




vector<double> getMagnAscii(const double& _x,const double& _y,const double& _z){
    //Interpolate potentials on r,z direction and calculate stuff
     //Currently works only if r>0 and z>0.
     r=sqrt(_x*_x+_y*_y);if(r == 0)r=1e-10; //otherwise division by zero!     
     //check of ze wel in of bound zitten he!
     if(_x>Br_max||_y>Br_max||r>Br_max||_z>Bz_max||_z<(Bz_min+Bz_step)){
           cout<<"B particle out of bound: "<<_x<<" "<<_y<<" "<<_z<<endl;}
     //scale both particles
      z_index= _z-Bz_min; 
      r_index= r-Br_min;

     
      r1_index=fabs(floor(r_index/Br_step));r2_index=r1_index+1;//ceil(r_index/Br_step); //the next one he: (r+0.5*r_step)/r_step
      z1_index=floor(z_index/Bz_step);z2_index=z1_index+1;//ceil(z_index/Bz_step);
    //  cout<<r1_index<<"\t"<<r2_index<<endl;
    //  cout<<z1_index<<"\t"<<z2_index<<endl;
      
     B12=(*Bmatrix)(z2_index,r1_index);B22=(*Bmatrix)(z2_index,r2_index);
     B11=(*Bmatrix)(z1_index,r1_index);B21=(*Bmatrix)(z1_index,r2_index);
      
      r1=r1_index*Br_step+Br_min;r2=r2_index*Br_step+Br_min;
      z1=z1_index*Bz_step+Bz_min;z2=z2_index*Bz_step+Bz_min;
      
   
      //interpolate Potential!  http://en.wikipedia.org/wiki/Bilinear_interpolation
      dzr=1.0/((r2-r1)*(z2-z1));
      Bzr.first = B11.first*(r2-r)*(z2-_z)*dzr
                 +B21.first*(r-r1)*(z2-_z)*dzr
                 +B12.first*(r2-r)*(_z-z1)*dzr
                 +B22.first*(r-r1)*(_z-z1)*dzr;
      Bzr.second = B11.second*(r2-r)*(z2-_z)*dzr
                 +B21.second*(r-r1)*(z2-_z)*dzr
                 +B12.second*(r2-r)*(_z-z1)*dzr
                 +B22.second*(r-r1)*(_z-z1)*dzr;
 	//cout<<"asked for B with z="<<_z<<" r="<<r<<"  ->Bz="<<Bzr.first<<" // Br="<<Bzr.second<<" T"<<endl;cin.get();
      if(Bzr.second != Bzr.second){ 
          //  cout<<"Bzr.second "<<Bzr.second<<endl;
            cout<<"Bfoutje: "<<dzr<<endl;
            cout<<"x: "<<_x<<endl;
            cout<<"y: "<<_y<<endl;
            cout<<"r: "<<r<<endl;
            cout<<"z: "<<_z<<endl; 
      cout<<r1_index<<"\t"<<r2_index<<endl;
      cout<<z1_index<<"\t"<<z2_index<<endl;               
            cin.get();
            Bzr.second = 0.0;
            magnetic_field[0]=0.0; //B on position 0.0 he!
            magnetic_field[1]=0.0;
      }
      else{
          /*if (fabs(_x) < 0.0001)*/ magnetic_field[0]=Bzr.second*(_x/r); //else magnetic_field[0]=0.0;
         /* if (fabs(_y) < 0.0001)*/ magnetic_field[1]=Bzr.second*(_y/r); //else magnetic_field[1]=0.0;
     }
     
     if(Bzr.first != Bzr.first){ 
         //  cout<<"Bzr.first "<<Bzr.first<<endl;        
           Bzr.first=0.0;
           magnetic_field[2]=0.0;                                  
     }
     else{
          magnetic_field[2]=Bzr.first;
     }
  
     return magnetic_field;                
}




/*    
void ReadVoltages(char *_V_filename){
        
     ifstream V_str;
     V_str.open(_V_filename); 
     if (!V_str){cout<<"Error in opening "<<_V_filename<<". The program will be stopped."<<endl; exit(1);}
     
     double r,z,V;
     int r_stepper = 0, z_stepper = 0;

     char firstline[256]; V_str.getline(firstline, 256);
     //read file minus first whiteline.
     bool valuesset = false;
     while ( !V_str.eof() && z_stepper != Vmatrix.nrows() && r_stepper != Vmatrix.ncols() ) { // keep reading until end-of-file
        V_str>>r>>z>>V;
        //r_min,max en voor z ook bepalen he.
        if(!valuesset){valuesset=true;r_min=r_max=r;z_min=z_max=z;}
        if(r<r_min){r_min=r;}
        if(r>r_max){r_max=r;}
        if(z<z_min){z_min=z;}
        if(z>z_max){z_max=z;}        
        Vmatrix(z_stepper,r_stepper)=V;
        z_stepper++;
        if (z_stepper == Vmatrix.nrows()){z_stepper =0;r_stepper++;}                      
     }
     //set the r and z_step
     r_step=(r_max-r_min)/(Vmatrix.ncols()-1);
     z_step=(z_max-z_min)/(Vmatrix.nrows()-1);
    
    //  cout<<r_min<<" "<<r_max<<" "<<r_step<<endl;
    //  cout<<z_min<<" "<<z_max<<" "<<z_step<<endl;
     
     cout<<"file: "<<_V_filename<<" read\n";
    // cout<<"*-*"<<Vmatrix(55,55)<<"*-*\n";
     V_str.close(); 
} 
*/
  
  
  
/*
vector<double>& getElectr_via_voltage(const double& _x,const double& _y,const double& _z){
     double x,y,z;
     x = _x;y = _y;z= _z;
     r=sqrt(x*x+y*y);if(r == 0)r=1e-10; //otherwise division by zero!     
   
     
     //check of ze wel in of bound zitten he!
  //   if(x>r_max||y>r_max||z>z_max||x<r_min||y<r_min||z<z_min){
  //         cout<<"particle out of bound: "<<x<<" "<<y<<" "<<z<<endl;}
     //scale both particles he.
     //uppercorner is(z_min,r_min)
      z= z+(-z_min); 
      r= r+(-r_min);

     
     //find correct neighbours in array
      double electric_fieldr, electric_fieldz;
      double rs[4],zs[4];
      double volts[4];
      int coord[4];
      
      rs[0] = fabs(r - r_step*0.5);
      rs[1] = r + r_step*0.5;
      rs[2] = rs[3] = r;
      
      zs[2] = z + z_step*0.5;
      zs[3] = z - z_step*0.5;
      zs[0] = zs[1] = z;
      
      double u = 0, t = 0;
      
      for (int i = 0;i<4;i++){
          coord[0] = (int)(rs[i]/r_step);
          coord[1] = (int)(rs[i]/r_step)+1;
          coord[2] = (int)(zs[i]/z_step);
          coord[3] = (int)(zs[i]/z_step)+1;
      
          t = (rs[i]-coord[0]*r_step)/r_step;
          u = (zs[i]-coord[2]*z_step)/z_step;
          
         if((coord[1] >= Vmatrix.ncols())|| (coord[3]>=Vmatrix.nrows()) || coord[2]<0.0 || coord[0]<0.0){
                cout<<"Array at of bound error in 'void getElectr(double x,double y,double z)'..."<<endl;
                cout<<"rs[i] = "<<rs[i]<<endl;
                cout<<"coord[1] = "<<coord[1]<<endl;
                cout<<"zs[i] = "<<zs[i]<<endl;
                cout<<"coord[3] = "<<coord[3]<<endl;
                cout<<"Program terminated"<<endl;   
                exit(1); 
                break;  
          }
          volts[i] = (1-t)*(1-u)*Vmatrix(coord[2],coord[0])+
                t*(1-u)*Vmatrix(coord[2],coord[1])+
                t*u*Vmatrix(coord[3],coord[1])+
                (1-t)*u*Vmatrix(coord[3],coord[0]);
      }
      electric_fieldr = (volts[0]-volts[1])/r_step;
      electric_fieldz = (volts[3]-volts[2])/z_step;


      
      
      electric_field[0] = electric_fieldr*(x/r);
      electric_field[1] = electric_fieldr*(y/r);
      electric_field[2] = electric_fieldz;
      return electric_field;
 
}
*/


/*
vector<double>& getElectrInterpolation(const double& _x,const double& _y,const double& _z){
     //Interpolate potentials on r,z direction and calculate stuff
     //Currently works only if r>0 and z>0.
     
     double x,y,z,r;
     x = _x;y = _y;z= _z;
     r=sqrt(x*x+y*y);if(r == 0)r=1e-10; //otherwise division by zero!     
   
     
     //check of ze wel in of bound zitten he!
     if(x>r_max||y>r_max||z>z_max||z<z_min){
           cout<<"particle out of bound: "<<x<<" "<<y<<" "<<z<<endl;cin.get();}
     //scale both particles
      double z_index= z+(-z_min); 
      double r_index= r+(-r_min);

      double z1,z2,r1,r2;
      double z1_index,z2_index,r1_index,r2_index;

      r1_index=(int)(r_index/r_step-0.5);r2_index=(int)(r_index/r_step+0.5); //the next one he: (r+0.5*r_step)/r_step
      z1_index=(int)(z_index/z_step-0.5);z2_index=(int)(z_index/z_step+0.5);
      
      //matrix is first row, then column he.
      double V12=Vmatrix(z2_index,r1_index);double V22=Vmatrix(z2_index,r2_index);
      double V11=Vmatrix(z1_index,r1_index);double V21=Vmatrix(z1_index,r2_index);
      
      r1=r1_index*r_step+r_min;r2=r2_index*r_step+r_min;
      z1=z1_index*z_step+z_min;z2=z2_index*z_step+z_min;
      
   
      //interpolate Potential!  http://en.wikipedia.org/wiki/Bilinear_interpolation
      double dz=z2-z1;
      double dr=r2-r1;
      double dzr=dr*dz;
      double Vzr= V11*(r2-r)*(z2-z)/dzr
                 +V21*(r-r1)*(z2-z)/dzr
                 +V12*(r2-r)*(z-z1)/dzr
                 +V22*(r-r1)*(z-z1)/dzr;

      cout.precision(8);*/
      /*   r1,z2----------Vup----r2,z2
            |              |      |      
            |              |      |
          Vleft-----------Vrz----Vright                        
            |              |      |
            |              |      |
          r1,z1-----------Vdown--r2,z1        */
/*      double Vleft,Vright,Vup,Vdown;
      //interpolate potentials  
      Vleft =V11+(z-z1)*(V12-V11)/(z2-z1);
      Vright=V21+(z-z1)*(V22-V21)/(z2-z1);      
      Vup   =V12+(r-r1)*(V22-V12)/(r2-r1);
      Vdown =V11+(r-r1)*(V21-V11)/(r2-r1);
     //Calculate Electric Field left and right
      double Erleft=-(Vzr-Vleft)/(r-r1);
      double Erright=-(Vright-Vzr)/(r2-r);
      double Ezup=-(Vzr-Vup)/(z-z2);
      double Ezdown=-(Vdown-Vzr)/(z1-z);
      //Interpolate Electric Field
      double Er=Erleft+(r-r1)*(Erright-Erleft)/(r2-r1);
      double Ez=Ezup+(z-z2)*(Ezdown-Ezup)/(z1-z2);

      
 */     //methode a2) of via interpolatie Voltages op de zijassen berekenen, en dan enkel dichtste bij punten nemen!
   /*   Er=(Vright-Vzr)/(r2-r);
      Ez=(Vup-Vzr)/(z2-z);
   */   
      //methode b) of gewoon meer gemiddeldes nemen op de hoekpunten van omringen vierkand
   /*   
    */  
      //methode c) gewoon via hoekpunten, redelijk efficient!!
  /*    double Er=(-(V22-V12)/r_step-(V21-V11)/r_step)*0.5;
      double Ez=(-(V12-V11)/z_step-(V22-V21)/z_step)*0.5;
      double V=Vzr;   
   */ 
      
 /*     electric_field[0]=Er*(x/r);
      electric_field[1]=Er*(y/r);
      electric_field[2]=Ez;
      return electric_field;
}*/

/*
vector<double>& getElectrFast(const double& _x,const double& _y,const double& _z){
     //Interpolate potential+get potential on closest point.
     //Calculate Electric field due to this 2 potentials.
     //Not yet implemented
     
     double x,y,z,r;
     x = _x;y = _y;z= _z;
     r=sqrt(x*x+y*y);if(r == 0)r=1e-10; //otherwise division by zero!     
   
     
     //check of ze wel in of bound zitten he!
     if(x>r_max||y>r_max||z>z_max||z<z_min){
           cout<<"particle out of bound: "<<x<<" "<<y<<" "<<z<<endl;cin.get();}
     //scale both particles
      double z_index= z+(-z_min); 
      double r_index= r+(-r_min);

      double z1,z2,r1,r2;
      double z1_index,z2_index,r1_index,r2_index;

      //to be implemented
      
  //    electric_field[0]=Er*(x/r);
  //    electric_field[1]=Er*(y/r);
  //    electric_field[2]=Ez;
      return electric_field;
}

*/


/*

      Very old function

*/

/*
void ConvertBfile(char *_B_filename){
    ofstream testie;
    testie.open("9T01T_converted.txt"); 

    int x_scale_mags = 49, z_scale_mags = 1161;
    double grid_scale_mag = 0.0025;//0.0001 decay cooler 0.00005 doubled //m/gu
    double grid_shift_magn = -1.106; //in m, compared to ions coordinates

    
     double x,z,B;
     ifstream magneticfield;
     magneticfield.open(_B_filename); 
     if (!magneticfield)
     {cout<<"Error in opening magnet file. The program will be stopped."<<endl; exit(1);}
     
     for (int i = 0;i<z_scale_mags;i++){
         for (int j = 0;j<x_scale_mags;j++){
              magneticfield>>B;
              z=(1106-i)*0.25;//veelkans nog schalen en verschuiven
              x=20.0/50.0*j;//20.0/50.0=0.04 dus ongeveer 0.25^-1
              testie<<x<<"\t"<<z<<"\t"<<B<<endl;
         }         
     }
    testie.close(); 
    magneticfield.close();
    
    cout<<"Finished reading MagneticField"<<endl;   
     
}
   

     
void getNearestNeighbours(const double& _x,const double& _y,const double& _z){
     double x,y,z,r;
     x = _x;y = _y;z= _z;
        //check of ze wel in of bound zitten he!
     if(x>r_max||y>r_max||z>z_max||x<r_min||y<r_min||z<z_min){
           cout<<"particle out of bound: "<<x<<" "<<y<<" "<<z<<endl;}

      int rs[4],zs[4];
      double r_scale=1.0/r_step;
      double z_scale=1.0/z_step;
      

      r = sqrt(x*x + y*y);
      if(r == 0)r=1e-10;

//scale both particles he.
//uppercorner is(z_min,r_min)
      r = r+(-r_min);
      z= z+(-z_min);
      
           rs[0]=(int)(r/r_step);
           rs[1]=(int)(r/r_step-0.5); //the next one he: (r+0.5*r_step)/r_step
           rs[2]=(int)(r/r_step+0.5);
           zs[0]=(int)(z/z_step);
           zs[1]=(int)(z/z_step-0.5);
           zs[2]=(int)(z/z_step+0.5);
      cout<<"************"<<endl;
      //cout<<Vmatrix(rs[0],z);    
      cout<<"Voltage@ "<<sqrt(_x*_x+_y*_y)<<","<<_z<<" is: "<<Vmatrix(zs[0],rs[0])<<"\n";
    //  cout<< Vmatrix(zs[1],rs[1])<<" "<<Vmatrix(zs[1],rs[2])<<endl;
    //  cout<< Vmatrix(zs[2],rs[1])<<" "<<Vmatrix(zs[2],rs[2])<<endl;   
}     
  */   
