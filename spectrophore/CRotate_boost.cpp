#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <boost/python.hpp>

using namespace std; 

// Prototypes
//exposed funcitions
boost::python::object CRotate(boost::python::object rotationStepList, boost::python::object coor, boost::python::object prop, 
								std::string normalization, float resolution, int beginProbe, int endProbe, boost::python::object radii) ; 

//hidden functions
vector<double> Calculate(vector< vector<double> >  Oricoor, int size_coor, vector< vector<double> > Property,
							string normalization, vector<int> rotationStepList, float resolution, int beginProbe, int endProbe, vector<double> radius ); 

vector< vector<double> > getBox (vector< vector<double> > Matrix1, int size_coor, float resolution, vector<double> radius ); 

vector<double> getEnergy (vector< vector<double> > boxPoint, vector< vector<double> > coor, int size_coor, vector< vector<double> > property,
							 int beginProbe, int endProbe);

vector< vector<double> > dotMatrix (vector< vector<double> > oricoor, const int size_coor, const double rm[3][3]); // rm = rotation matrix

// End of Protorypes
boost::python::object CRotate(boost::python::object rotationStepList, boost::python::object coor, 
	boost::python::object prop, std::string normalization, float resolution, int beginProbe, int endProbe, boost::python::object radii) {
	//setting variables
	int _size = len(rotationStepList);
	int _size_coor = len(coor);

	// Declaring list
	vector<int> v(_size);

	for (int i = 0 ; i < _size ; i++) {
		v[i] = boost::python::extract<int>(rotationStepList[i]);
	}	

	vector< vector<double> > points;
	vector< vector<double> > property;

	for (int i = 0; i < _size_coor; i++) {
		boost::python::object rowlist = boost::python::extract<boost::python::object>(coor[i]);

		vector<double> row;
		for (int j = 0; j < 3; j++) {
			row.push_back(boost::python::extract<double>(rowlist[j]));
		}
		points.push_back(row);

		//////////////////
		boost::python::object rowProp = boost::python::extract<boost::python::object>(prop[i]);

		vector<double> row1;
		for (int j = 0; j < 4; j++) {	
			row1.push_back(boost::python::extract<double>(rowProp[j]));
		}
		property.push_back(row1);		
	}

	// getting radii 
	vector<double> radius(_size_coor);

	for (int i = 0 ; i < _size_coor ; i++) {
		radius[i] = boost::python::extract<double>(radii[i]);
	}

	vector<double> spectrophore =  Calculate(points, _size_coor, property, normalization, v, resolution, beginProbe, endProbe, radius); // vector 2D, int vector2D, vector

	boost::python::list l;

	for (unsigned int i = 0; i < spectrophore.size(); i++) {
		l.append(spectrophore[i]);
	}

	return l;
}

vector< vector<double> > getBox (vector< vector<double> > Matrix1, int size_coor, float resolution, vector<double> radius ) {
	double m[3] = {};
	double p[3] = {};
	double h[3] = {};

	for (int i = 0; i < size_coor; i++) {
		for (int j = 0; j < 3; j++) {
			if ((Matrix1[i][j] - resolution - radius[i]) < m[j]) 
				m[j] = Matrix1[i][j] - resolution - radius[i];

			if ((Matrix1[i][j] + resolution + radius[i]) > p[j])
				p[j] = Matrix1[i][j] + resolution + radius[i];
		}
	}

	for (int i = 0; i < 3; ++i) {
		h[i] = (m[i] + p[i])/2;
	}

	double temp_boxpoint[3][12] = {{h[0],p[0],h[0],m[0],m[0],p[0],m[0],p[0],p[0],h[0],m[0],h[0]},
						           {m[1],h[1],p[1],h[1],m[1],m[1],p[1],p[1],h[1],m[1],h[1],p[1]},
						           {p[2],p[2],p[2],p[2],h[2],h[2],h[2],h[2],m[2],m[2],m[2],m[2]}};


	vector< vector<double> > boxpoint(12, vector<double>(7, 0));
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 3; j++) {
			boxpoint[i][j] = temp_boxpoint[j][i];
		}
	}

	return boxpoint;
}

vector<double> Calculate (vector< vector<double> >  Oricoor, int size_coor, vector< vector<double> > Property, string normalization, 
	vector<int> rotationStepList, float resolution, int beginProbe, int endProbe, vector<double> radius) {
	const int N_PROPERTIES = 4;	
	int size = rotationStepList.size();	
	const int numberOfProbes = endProbe - beginProbe;

    //self.__getBox__(oricoor)
	vector< vector<double> > boxPoint = getBox(Oricoor, size_coor, resolution, radius);

    //self.__getEnergies__(oricoor, prop)
	vector<double> e = getEnergy(boxPoint, Oricoor, size_coor, Property, beginProbe, endProbe);

    //self.sphore = self.energy
	vector<double> sphore(N_PROPERTIES * numberOfProbes);
    for (int i = 0; i < N_PROPERTIES * numberOfProbes; ++i) {
        sphore[i] = e[i];
    }

	for (int i = 0; i < size; i++) {
		int rotationStep = rotationStepList[i];
		//cout << "Rotation Sterp: " << rotationStep << endl;
		for (int iTheta = 0; iTheta < 180; iTheta += rotationStep){
			double theta = 0.017453292519943 * iTheta;
			double cos_theta = cos(theta);
			double sin_theta = sin(theta);

			double rotateY[3][3] = {{cos_theta, 0, -sin_theta},
									{        0, 1,          0},
									{sin_theta, 0, cos_theta}};

			vector< vector<double> > Matrix_Y = dotMatrix(Oricoor, size_coor, rotateY);

			for (int iPsi = 0; iPsi < 360; iPsi += rotationStep) {
				double psi = 0.017453292519943 * iPsi;
			    double cos_psi = cos(psi);
			    double sin_psi = sin(psi);

				double rotateZ[3][3] = {{cos_psi, -sin_psi, 0},
										{sin_psi,  cos_psi, 0},
										{      0,        0, 1}};

				vector< vector<double> > Matrix_Z = dotMatrix(Matrix_Y, size_coor, rotateZ);

	    		for (int iPhi = 0; iPhi < 360; iPhi += rotationStep){
					double phi = 0.017453292519943 * iPhi;
					double cos_phi = cos(phi);
					double sin_phi = sin(phi); 

					double rotateX[3][3] = {{1,       0,        0},
											{0, cos_phi, -sin_phi},
											{0, sin_phi,  cos_phi}};

					vector< vector<double> > Matrix_X = dotMatrix(Matrix_Z, size_coor, rotateX);
					vector< vector<double> > boxP  = getBox(Matrix_X, size_coor, resolution, radius);

				    //self.__getEnergies__(oricoor, prop)
					vector<double> energy = getEnergy(boxP, Matrix_X, size_coor, Property, beginProbe, endProbe); 

				    //self.sphore = self.energy
				    for (int a = 0; a < N_PROPERTIES * numberOfProbes; ++a) {
				    	if (energy[a] < sphore[a]) {
				        	sphore[a] = energy[a];
				        }
				    }
				}
			}
		}
	}

    vector<double> spectrophore(N_PROPERTIES * numberOfProbes);

    for (int i = 0; i < N_PROPERTIES * numberOfProbes; i++) {
    	spectrophore[i] = - (100 * sphore[i]);
    }

   // Normalisation
   double mean[N_PROPERTIES];
   double std[N_PROPERTIES];
   unsigned int m;
   for (int i(0); i < N_PROPERTIES; ++i) {
      mean[i] = 0.0;
      for (int n(0); n < numberOfProbes; ++n) {
         m = (i * numberOfProbes) + n;
         mean[i] += spectrophore[m];
      }


      mean[i] /= float(numberOfProbes);
      std[i] = 0.0;
      for (int n(0); n < numberOfProbes; ++n) {
         m = (i * numberOfProbes) + n;
         std[i] += (spectrophore[m] - mean[i]) * (spectrophore[m] - mean[i]);
      }
      std[i] /= float(numberOfProbes);
      std[i] = sqrt(std[i]);
   }

   if ((normalization == "mean") || (normalization == "all")) {
      for (int i(0); i < N_PROPERTIES; ++i) {
         for (int n(0); n < numberOfProbes; ++n) {
            m = (i * numberOfProbes) + n;
            spectrophore[m] -= mean[i];
         }
      }
   }
   if ((normalization == "std") || (normalization == "all")) {
      for (int i(0); i < N_PROPERTIES; ++i) {
         for (int n(0); n < numberOfProbes; ++n) {
            m = (i * numberOfProbes) + n;
            spectrophore[m] /= std[i];
         }
      }
   }

   // Return
   return spectrophore; 
}

vector< vector<double> > dotMatrix (vector< vector<double> > oricoor, const int size_coor, const double rm[3][3]){
	vector< vector<double> > newMatrix(size_coor , vector<double>(3, 0));
	//double newMatrix[size_coor][3];
	for (int i = 0; i < size_coor; i++)	{
		for (int j = 0; j < 3; j++) {
			newMatrix[i][j] =   oricoor[i][0]*rm[0][j] +
								oricoor[i][1]*rm[1][j] +
								oricoor[i][2]*rm[2][j] ;
		}
	}

	return newMatrix;
}

vector<double> getEnergy (vector< vector<double> > boxPoint, vector< vector<double> > coor, int size_coor, vector< vector<double> > property, int beginProbe, int endProbe) {
	const int N_PROPERTIES = 4;
	const int numberOfProbes = endProbe - beginProbe;

	int _probe[48][12] = {{1, 1,-1,-1,-1, 1, 1,-1,-1,-1, 1, 1},
			 {1, 1,-1,-1, 1,-1,-1, 1,-1,-1, 1, 1},
			 {1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1, 1},
			 {1, 1, 1,-1,-1,-1,-1,-1, 1, 1,-1, 1},
			 {1, 1, 1,-1,-1, 1,-1, 1,-1,-1, 1,-1},
			 {1, 1, 1,-1, 1,-1, 1,-1,-1,-1, 1,-1},
			 {1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1,-1},
			 {1, 1, 1, 1,-1,-1,-1,-1, 1,-1, 1,-1},
			 {1, 1, 1, 1,-1,-1,-1,-1, 1, 1,-1,-1},
			 {1, 1, 1, 1, 1,-1,-1, 1,-1,-1,-1,-1},
			 {1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1},
			 {1, 1, 1,-1,-1, 1,-1,-1,-1, 1,-1, 1},
			 {1, 1,-1,-1,-1,-1, 1, 1,-1, 1, 1,-1},
			 {1, 1, 1,-1,-1,-1,-1,-1, 1,-1, 1, 1},
			 {1, 1, 1,-1,-1,-1,-1, 1,-1, 1, 1,-1},
			 {1, 1, 1,-1,-1,-1, 1,-1,-1, 1,-1, 1},
			 {1, 1, 1,-1,-1,-1, 1,-1,-1, 1, 1,-1},
			 {1, 1, 1,-1,-1,-1, 1,-1, 1,-1, 1,-1},
			 {1, 1, 1,-1,-1,-1, 1,-1, 1, 1,-1,-1},
			 {1, 1, 1,-1,-1,-1, 1, 1,-1, 1,-1,-1},
			 {1, 1, 1,-1,-1,-1, 1, 1, 1,-1,-1,-1},
			 {1, 1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1},
			 {1, 1, 1,-1,-1, 1, 1,-1,-1,-1, 1,-1},
			 {1, 1, 1,-1,-1, 1, 1,-1,-1, 1,-1,-1},
			 {1, 1, 1,-1,-1, 1, 1, 1,-1,-1,-1,-1},
			 {1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1,-1},
			 {1, 1, 1,-1, 1,-1, 1,-1,-1,-1,-1, 1},
			 {1, 1, 1,-1, 1,-1, 1, 1,-1,-1,-1,-1},
			 {1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1},
			 {1, 1, 1, 1, 1,-1,-1,-1,-1,-1, 1,-1},
			 {1, 1,-1,-1, 1,-1, 1,-1, 1,-1,-1, 1},
			 {1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1},
			 {1, 1, 1,-1,-1, 1,-1,-1,-1,-1, 1, 1},
			 {1, 1, 1,-1, 1,-1,-1,-1,-1, 1,-1, 1},
			 {1, 1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1},
			 {1, 1, 1,-1, 1,-1,-1,-1, 1,-1, 1,-1},
			 {1, 1, 1,-1, 1,-1,-1,-1, 1,-1,-1, 1},
			 {1, 1, 1,-1, 1, 1,-1,-1,-1,-1,-1, 1},
			 {1, 1, 1,-1, 1, 1,-1,-1, 1,-1,-1,-1},
			 {1, 1, 1,-1, 1,-1,-1, 1,-1, 1,-1,-1},
			 {1, 1, 1,-1, 1,-1,-1, 1,-1,-1, 1,-1},
			 {1, 1, 1,-1, 1,-1,-1, 1,-1,-1,-1, 1},
			 {1, 1, 1,-1, 1, 1,-1, 1,-1,-1,-1,-1},
			 {1, 1, 1,-1, 1, 1,-1,-1,-1,-1, 1,-1},
			 {1, 1, 1,-1, 1,-1, 1,-1,-1, 1,-1,-1},
			 {1, 1, 1,-1, 1, 1, 1,-1,-1,-1,-1,-1},
			 {1, 1, 1, 1, 1,-1,-1,-1, 1,-1,-1,-1},
			 {1, 1, 1, 1, 1,-1,-1,-1,-1, 1,-1,-1}};

	double distance;

	vector<double> e(N_PROPERTIES * numberOfProbes);
	//double e[N_PROPERTIES * numberOfProbes] = {};
   // Distance between atom and each boxpoint
	for (unsigned int point = 0; point < 12; point++) {
      // Reset value at box point
		for (unsigned int prop = 3; prop < 7; prop++) {
			boxPoint[point][prop] = 0.0;
		}

		for (int i = 0; i < size_coor; i++) {
			distance = sqrt(pow((coor[i][0] - boxPoint[point][0]), 2) +
								  pow((coor[i][1] - boxPoint[point][1]), 2) +
								  pow((coor[i][2] - boxPoint[point][2]), 2));

			for (unsigned int prop = 0; prop < N_PROPERTIES; prop++) {
				boxPoint[point][prop + 3] += property[i][prop]/distance; /// tem que passar as propriedades para os i atomos
			}			
		}
	}

   // Reset energy to zero
   for ( int i = 0; i < N_PROPERTIES * numberOfProbes; i++) {
      e[i] = 0.0;
   }

   unsigned int index;
   for (unsigned int prop = 0; prop < N_PROPERTIES; prop++) {
      for (unsigned int point = 0; point < 12; point++) {
         // Calculate for each probe the total interaction energy
         for ( int probe = beginProbe; probe < endProbe; probe++) {
            index = prop * numberOfProbes + (probe - beginProbe);
            e[index] += boxPoint[point][prop +3] *_probe[probe][point];
         }
      }
   }

    return e;
}

BOOST_PYTHON_MODULE(CRotate){
	def("CRotate", &CRotate);
}
