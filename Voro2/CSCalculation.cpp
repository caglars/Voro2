//
//  CSCalculation.cpp
//  Voro2
//
//  Created by Caglar Sinayuc on 12/09/14.
//  Copyright (c) 2014 caglars. All rights reserved.
//
#include <ctime>
#include "CSCalculation.h"
#include "voro++.hh"
#include "CSSolver.h"
#include <math.h>
using namespace voro;

#include <vector>
using namespace std;

// Set up constants for the container geometry
const double x_min=0,x_max=6000;
const double y_min=0,y_max=6000;
const double z_min=0,z_max=120;
const double cvol=(x_max-x_min)*(y_max-y_min)*(z_max-z_min);

// Set up the number of blocks that the container is divided into
const int n_x=60,n_y=60,n_z=6;

// The sampling distance for the grids of find_voronoi_cell calls
//const double h=0.05;

// The cube of the sampling distance, corresponding the amount of volume
// associated with a sample point
//const double hcube=h*h*h;

// Set the number of particles that are going to be randomly introduced
const int particles=18;

// This function returns a random double between 0 and 1
double rnd() {
    return double(rand())/RAND_MAX;
}

//container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);


void Calculation::calculate() {
    int i;
	double x,y,z,r;
	vector<int> neigh, f_vert;
	vector<double> vec, faceAreas;
    srand( (unsigned)time( NULL ) );
    
    unsigned long maxNumberOfNeighbors = 0; // henuz bir yerde kullanmadim belki lazim olur
    
    container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);
    
/*
    for(i=0;i<particles-1;i++) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=40;
		con.put(i,x,y,z);
	}
    con.put(particles-1, 3000, 3000, 40);
*/
    

    con.put(0, 1000, 1000, 40);
    con.put(1, 3000, 1000, 40);
    con.put(2, 5000, 1000, 40);
    con.put(3, 1000, 3000, 40);
    con.put(4, 5000, 5000, 40);
    con.put(5, 5000, 3000, 40);
    con.put(6, 1000, 5000, 40);
    con.put(7, 3000, 5000, 40);
    con.put(8, 3000, 3000, 40);
 
    con.put(9, 1000, 1000, 80);
    con.put(10, 3000, 1000, 80);
    con.put(11, 5000, 1000, 80);
    con.put(12, 1000, 3000, 80);
    con.put(13, 5000, 5000, 80);
    con.put(14, 5000, 3000, 80);
    con.put(15, 1000, 5000, 80);
    con.put(16, 3000, 5000, 80);
    con.put(17, 3000, 3000, 80);

    
    
    /*
    con.put(0, 550, 1000, 40);
    con.put(1, 2300, 1200, 40);
    con.put(2, 5400, 1100, 40);
    con.put(3, 1200, 3500, 40);
    con.put(4, 2000, 3000, 40);
    con.put(5, 6300, 3200, 40);
    con.put(6, 1400, 5500, 40);
    con.put(7, 2500, 6000, 40);
    con.put(8, 6000, 6000, 40);
     */
    
    //int particles = 0;
    
    positionX = new double[particles];
    positionY = new double[particles];
    positionZ = new double[particles];
    
    //double positionX[particles];
    //double positionY[particles];
    //double positionZ[particles];
    for (int count = 0; count < particles; count++) {
        positionX[count]=0;
        positionY[count]=0;
        positionZ[count]=0;
    }
    
    //double cellVolume[particles];
    cellVolume = new double[particles];
    for (int counter=0; counter<particles; counter++) {
        cellVolume[counter] = 0;
    }
    
    // Sum up the volumes, and check that this matches the container volume
	double vvol=con.sum_cell_volumes();
	printf("Container volume : %g\n"
	       "Voronoi volume   : %g\n"
	       "Difference       : %g\n",cvol,vvol,vvol-cvol);
	
    FILE *f1=safe_fopen("find_voro_cell.vol","w");
	FILE *f2=safe_fopen("find_voro_cell_v.gnu","w");
	c_loop_all cla(con);
	voronoicell_neighbor v;
	if(cla.start()) do if (con.compute_cell(v,cla)) {
    
        
		// Get the position and ID information for the particle
		// currently being considered by the loop. Ignore the radius
		// information.
		cla.pos(i,x,y,z,r);
        
        positionX[i] = x;
        positionY[i] = y;
        positionZ[i] = z;
        cellVolume[i] = v.volume();
        
        
		v.neighbors(neigh);
		//printf("particle: %d, neighbors: %lu\n", i, neigh.size());
		//printf("neighbors: %d, %d, %d, %d, %d, %d\n", neigh[0], neigh[1], neigh[2], neigh[3], neigh[4], neigh[5]);
        
		// For each face, this routine outputs a bracketed sequence of numbers
		// containing a list of all the vertices that make up that face.
		v.face_vertices(f_vert);
		//printf("particle: %d, faces: %lu\n", i, f_vert.size());
		//printf("faces: %d, %d, %d, %d, %d, %d\n", f_vert[0], f_vert[1], f_vert[2], f_vert[3], f_vert[4], f_vert[5]);
        
		v.vertices(x,y,z,vec);
		//printf("particle: %d, x: %f, y: %f, z: %f, vertices: %lu\n", i, x, y, z, v.size());
		//printf("vertices: %f, %f, %f, %f, %f, %f\n", v[0], v[1], v[2], v[3], v[4], v[5]);
		//printf("number of vertices %d\n", c.p);
        
        //v.face_areas(faceAreas);
        
        if (neigh.size() > maxNumberOfNeighbors) maxNumberOfNeighbors = neigh.size();
        
        
        
        
		// Output the Voronoi cell to a file in gnuplot format
        v.draw_gnuplot(0,0,0,"simple_cell.gnu");
        
        printf("Cell ID             :%d\n", cla.pid());
        
        // Output vertex-based statistics
        printf("Total vertices      : %d\n",v.p);
        printf("Vertex positions    : ");v.output_vertices();puts("");
        printf("Vertex orders       : ");v.output_vertex_orders();puts("");
        printf("Max rad. sq. vertex : %g\n\n",0.25*v.max_radius_squared());
        
        // Output edge-based statistics
        printf("Total edges         : %d\n",v.number_of_edges());
        printf("Total edge distance : %g\n",v.total_edge_distance());
        printf("Face perimeters     : ");v.output_face_perimeters();puts("\n");
        
        // Output face-based statistics
        printf("Total faces         : %d\n",v.number_of_faces());
        printf("Surface area        : %g\n",v.surface_area());
        printf("Face freq. table    : ");v.output_face_freq_table();puts("");
        printf("Face orders         : ");v.output_face_orders();puts("");
        printf("Face areas          : ");v.output_face_areas();puts("");
        printf("Face normals        : ");v.output_normals();puts("");
        printf("Face vertices       : ");v.output_face_vertices();puts("");
        printf("Neighbors           : ");v.output_neighbors();puts("\n");
        
        // Output volume-based statistics
        v.centroid(x,y,z);
        printf("Volume              : %g\n"
               "Centroid vector     : (%g,%g,%g)\n",v.volume(),x,y,z);
        
        

        
		// Draw the Voronoi cell
		v.draw_gnuplot(x,y,z,f2);
	} while (cla.inc());
	fclose(f1);
	fclose(f2);
    
	// Do a custom output routine to store the number of vertices, edges,
	// and faces of each Voronoi cell
	con.print_custom(
                     "ID=%i, pos=(%x,%y,%z), vertices=%w, edges=%g, faces=%s",
                     "random.custom1");
    
	// Do a custom output routine to store a variety of face-based
	// statistics. Store the particle ID and position, the number of faces
	// the total face area, the order of each face, the areas of each face,
	// the vertices making up each face, and the neighboring particle (or
	// wall) corresponding to each face.
	con.print_custom("%i %q %s %F %a %f %t %l %n","random.custom2");
    
	// Do a custom output routine that outputs the particle IDs and
	// positions, plus the volume and the centroid position relative to the
	// particle center
	con.print_custom("%i %q %v %c","random.custom3");
    
	// Also create POV-Ray output of the Voronoi cells for use in the
	// rendering
	con.draw_cells_pov("random_v.pov");
    
	// Output the particle positions in gnuplot format
	con.draw_particles("random_points_p.gnu");
    
	// Output the Voronoi cells in gnuplot format
	con.draw_cells_gnuplot("random_points_v.gnu");

    
    
    
    
    
    
    
    printf("number of particles: %d\n", particles);
    
    
    simRun(con);
    
    
    
    
    
    
    
    
    

    
}

double Calculation::distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt(pow((x2-x1),2) + pow((y2-y1),2) + pow((z2-z1),2));
}

void Calculation::simRun(container& aCon) {
    
    double permeability = 0.015;
    double viscosity = 10;
    double formationVolumeFactor = 1;
    double porosity = 0.18;
    double liquidCompressibility = 3.5E-6;
    double referansFormationVolumeFactor = 1;
    double alphaConstant = 5.615;
    double betaConstant = 1.127;
    double deltaTime = 1;
    double length = 0;
    double totalCoefficient = 0;
    
    int numberOfTimeSteps = 10;
    
    int i;
    double x,y,z,r;
	vector<int> neigh, f_vert;
	vector<double> vec, faceAreas;
    double productionRate[particles];
    for (int counter=0; counter<particles; counter++) {
        productionRate[counter] = 0;
    }
    // lets say we have a well producing at 150 stb/day from the last particle
    productionRate[particles-1] = -250;
    
    double gamma[particles];
    for (int counter=0; counter<particles; counter++) {
        gamma[counter] = 0;
    }
    
    double pressure[particles];
    //pressure = new double [particles];
    for (int counter=0; counter<particles; counter++) {
        pressure[counter]=6000;
    }
    
    
    
    double rightHandSide[particles];
    //rightHandSide = new double[particles];
    for (int counter=0; counter<particles; counter++) {
        rightHandSide[counter] = 0;
    }
    
   
    coefficient = new double *[particles];
    for(int i = 0; i <particles; i++)
        coefficient[i] = new double[particles];
    
    for (int counter1=0; counter1<particles; counter1++) {
        for (int counter2=0; counter2<particles; counter2++) {
            coefficient[counter1][counter2]=0.0;
        }
    }
    
    Solution mySimpleSolver;
    double *solutionArray;
    solutionArray = new double[particles];
    for (int counter; counter<particles; counter++) {
        solutionArray[counter] = 0.0;
    }
    
    FILE *f3=safe_fopen("matrixA.dat","w");
	FILE *f4=safe_fopen("matrixB.dat","w");
    FILE *f5=safe_fopen("matrixSolution.dat","w");
    FILE *f6=safe_fopen("datafile.dat","w");
    
    for(int timeSteps = 0; timeSteps < numberOfTimeSteps; timeSteps++) {

    
        c_loop_all myLoop(aCon);
        voronoicell_neighbor myCell;
        if(myLoop.start()) do if (aCon.compute_cell(myCell,myLoop)) {
            myLoop.pos(i,x,y,z,r);
            myCell.neighbors(neigh);
            myCell.face_areas(faceAreas);
            totalCoefficient = 0;
            
            
            printf("Particle: %d position: x: %f, y: %f, z: %f\n", i, x, y, z);
            printf("Particle: %d position: x: %f, y: %f, z: %f, volume: %f\n", i, positionX[i], positionY[i], positionZ[i], myCell.volume());
            for (int neighCounter = 0; neighCounter < neigh.size(); neighCounter++) {
                printf("Particle %d neighbor %d = %d face area = %f\n", myLoop.pid(), neighCounter, neigh[neighCounter], faceAreas[neighCounter]);
                
                if (neigh[neighCounter]>0) {
                    length = distance(positionX[i], positionY[i], positionZ[i], positionX[neigh[neighCounter]], positionY[neigh[neighCounter]], positionZ[neigh[neighCounter]]);
                    coefficient[i][neigh[neighCounter]]=betaConstant*faceAreas[neighCounter]*permeability/(viscosity*formationVolumeFactor*length);
                    totalCoefficient = totalCoefficient + coefficient[i][neigh[neighCounter]];
                    printf("Neigbor: %d position: x: %f, y: %f, distance: %f\n",neigh[neighCounter], positionX[neigh[neighCounter]], positionY[neigh[neighCounter]], length);
                }
                
                
            }
            
            gamma[i] = (cellVolume[i]*porosity*liquidCompressibility)/(alphaConstant*referansFormationVolumeFactor);
            coefficient[i][i] = -totalCoefficient-gamma[i]/deltaTime;
            
            rightHandSide[i] = -productionRate[i] - (gamma[i]/deltaTime)*pressure[i];
            printf("The i value is: %d and productionRate:%f gamma:%f pressure:%f rightHandSide: %f\n", i, productionRate[i], gamma[i], pressure[i], rightHandSide[i]);
            
            
        } while (myLoop.inc());
        
        
        
        for (int counter1=0; counter1<particles; counter1++) {
            //for (int counter2=0; counter2<particles; counter2++) {
            fprintf(f3, "%f %f %f %f %f %f %f %f %f\n",coefficient[counter1][0], coefficient[counter1][1], coefficient[counter1][2], coefficient[counter1][3], coefficient[counter1][4], coefficient[counter1][5], coefficient[counter1][6], coefficient[counter1][7], coefficient[counter1][8]);
            //}
        }
        
        
        for (int counter1=0; counter1<particles; counter1++) {
            fprintf(f4, "%f\n", rightHandSide[counter1]);
        }
        
        
        
        solutionArray = mySimpleSolver.simpleSolver(particles, coefficient, rightHandSide, pressure);
        
        for (int counter1=0; counter1<particles; counter1++) {
            pressure[counter1] = solutionArray[counter1];
        }
        
        fprintf(f5, "\nTime step: %d Day: %.f\n\n", timeSteps+1, (timeSteps+1)*deltaTime);
        fprintf(f6, "\nTime step: %d Day: %.f\n\n", timeSteps+1, (timeSteps+1)*deltaTime);
        
        for (int counter1=0; counter1<particles; counter1++) {
            fprintf(f5, "%f\n", solutionArray[counter1]);
        }
        
        for (int counter1=0; counter1<particles; counter1++) {
            fprintf(f6, "%f %f %f\n", positionX[counter1], positionY[counter1], solutionArray[counter1]);
        }
        
    } // for timeSteps
    
    fclose(f3);
	fclose(f4);
    fclose(f5);
	fclose(f6);
    
    
}