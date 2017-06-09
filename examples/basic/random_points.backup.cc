// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"


using namespace voro;

// Set up constants for the container geometry
const double x_min=-1,x_max=1;
const double y_min=-1,y_max=1;
const double z_min=-1,z_max=1;
const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=10;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i;
	double x,y,z;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Randomly add particles into the container
	for(i=0;i<particles;i++) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		con.put(i,x,y,z);
	}

	// Sum up the volumes, and check that this matches the container volume
	double vvol=con.sum_cell_volumes();

	printf("Container volume : %g\n"
	       "Voronoi volume   : %g\n"
	       "Difference       : %g\n",cvol,vvol,vvol-cvol);

	// Output the particle positions in gnuplot format
	//con.draw_particles("random_points_p.gnu");

	// Output the Voronoi cells in gnuplot format
	voronoicell vc;
	c_loop_all vl(con);
	//refer to interface example :)
	std::vector<int> v;
	std::vector<double> face_edges;
	if(vl.start()) do{
		con.compute_cell(vc,vl);
		vc.face_vertices(v);

	  printf("\r\ncomputing cell:\t");


  	//Amit: Refer to Common.cc line 62:
  // 	void voro_print_face_vertices(std::vector<int> &v,FILE *fp) {
  	//TODO:Validate this method
			int j,k=0,l;
			if(v.size()>0) {
				l=v[k++];
				if(l<=1) {
					if(l==1){
						printf("(%g,%g,%g)",vc.pts[v[k]-1]*0.5,vc.pts[v[k]]*0.5,vc.pts[v[k]+1]*0.5);

						k++;
					}
					else printf("()");
				} else {
					j=k+l;
					printf("(%g,%g,%g)",vc.pts[v[k]-1]*0.5,vc.pts[v[k]]*0.5,vc.pts[v[k]+1]*0.5); k++;
					while(k<j){
						printf("(%g,%g,%g)",vc.pts[v[k]-1]*0.5,vc.pts[v[k]]*0.5,vc.pts[v[k]+1]*0.5);
						k++;
					}
					printf(")");
				}
				while((unsigned int) k<v.size()) {
					l=v[k++];
					if(l<=1) {
						if(l==1){printf("(%g,%g,%g)",vc.pts[v[k]-1]*0.5,vc.pts[v[k]]*0.5,vc.pts[v[k]+1]*0.5); k++;}
						else printf(" ()");
					} else {
						j=k+l;
						printf("(%g,%g,%g)",vc.pts[v[k]-1]*0.5,vc.pts[v[k]]*0.5,vc.pts[v[k]+1]*0.5); k++;
						while(k<j){ printf("(%g,%g,%g)",vc.pts[v[k]-1]*0.5,vc.pts[v[k]]*0.5,vc.pts[v[k]+1]*0.5); k++;}
						printf(")");
					}
				}
			}
		//}
	}while(vl.inc());
	//Amit:now you can call vl.start() and vl.inc() in a do while loop

	printf("done");
	con.draw_cells_gnuplot("random_points_v.gnu");
	con.print_custom("order=%o, vertices=%p","vertices.test");

	con.print_custom(
                "ID=%i, pos=(%x,%y,%z), vertices=%w, edges=%g, faces=%s, face_verts=%t"
                ,"face.test");
}


