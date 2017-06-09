// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"


using namespace voro;

// Set up constants for the container geometry
const double x_min=-10,x_max=10;
const double y_min=-10,y_max=10;
const double z_min=-0.5, z_max=0.5;
const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=20;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

double* createPointArray(double x, double y , double z){
  double * point =  new double[3];
  point[0] = x;
  point[1] = y;
  point[2] = z;
  return point;
}

int main() {

  int * cellsPtr;
  double * verticesPtr;

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
    z=0;//+rnd()*(z_max-z_min);
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
  std::vector<double> vertices = std::vector<double>();
  std::vector<int> cells = std::vector<int>();

  std::vector<int> v;

  if(vl.start()) do{
    con.compute_cell(vc,vl);
    vc.face_vertices(v);
    //find the displacement vector
    double * pp = con.p[vl.ijk]+con.ps*vl.q;
    printf("\n//ALL VC POINTS:\n");
    printf("\nvc.Add(new List<Vector3>(){\n");
    double *ptsp=vc.pts;
    cells.push_back(vc.p *3);
    for(int i=0;i<3*vc.p;i+=3) {
      printf("new Vector3(%gf,%gf,%gf),\n",*pp+*(ptsp)*0.5,
      pp[1]+*(ptsp+1)*0.5,
      pp[2]+ *(ptsp+2)*0.5);
      vertices.push_back(*pp+*(ptsp++)*0.5);
      vertices.push_back(pp[1]+*(ptsp++)*0.5);
      vertices.push_back(pp[2]+ *(ptsp++)*0.5);
    }
    // for(int zz =0; zz < (vc.current_vertices-3); zz+=3){
    //   printf("new Vector3(%gf,%gf,%gf),\n", *pp+vc.pts[zz]*0.5,pp[1]+vc.pts[zz+1]*0.5,pp[2]+vc.pts[zz+2]*0.5 );
    // }
    printf("\n });\n");
    //printf("\nSIZE:%lu\n", v.size());

    /*
    //Amit: Refer to Common.cc line 62:
  //  void voro_print_face_vertices(std::vector<int> &v,FILE *fp) {
    //TODO:Validate this method
      int j,k=0,l,counter=0;
      if(v.size()>0) {

        //container.hh line 604
        //printf("%d",p[vl.ijk]+con.ps*vl.q);

        l=v[k++];
        if(l<=1) {
          if(l==1){
            //printf("new Vector3(%gf,%gf,%gf),\r\n",*pp+vc.pts[v[k]-1]*0.5,pp[1]+vc.pts[v[k]]*0.5,pp[2]+vc.pts[v[k]+1]*0.5);

            vertices.push_back(*pp+vc.pts[v[k]-1]*0.5);
            vertices.push_back(*pp+vc.pts[v[k]]*0.5);
            vertices.push_back(*pp+vc.pts[v[k]+1]*0.5);
            counter++;
            k++;
          }
          else{
            printf("()");
          }
        } else {
          j=k+l;
          //printf("new Vector3(%gf,%gf,%gf),\r\n",*pp+vc.pts[v[k]-1]*0.5,pp[1]+vc.pts[v[k]]*0.5,pp[2]+vc.pts[v[k]+1]*0.5);
          vertices.push_back(*pp+vc.pts[v[k]-1]*0.5);
          vertices.push_back(*pp+vc.pts[v[k]]*0.5);
          vertices.push_back(*pp+vc.pts[v[k]+1]*0.5);
          counter++;
          k++;

          while(k<j){
            //printf("new Vector3(%gf,%gf,%gf),\r\n",*pp+vc.pts[v[k]-1]*0.5,pp[1]+vc.pts[v[k]]*0.5,pp[2]+vc.pts[v[k]+1]*0.5);
            vertices.push_back(*pp+vc.pts[v[k]-1]*0.5);
            vertices.push_back(*pp+vc.pts[v[k]]*0.5);
            vertices.push_back(*pp+vc.pts[v[k]+1]*0.5);
            counter++;
            k++;
          }
          //printf(")");
        }

        while((unsigned int) k<v.size()) {

          l=v[k++];
          if(l<=1) {
            if(l==1){
              //printf("new Vector3(%gf,%gf,%gf),\r\n",*pp+vc.pts[v[k]-1]*0.5,pp[1]+vc.pts[v[k]]*0.5,pp[2]+vc.pts[v[k]+1]*0.5);
              vertices.push_back(*pp+vc.pts[v[k]-1]*0.5);
              vertices.push_back(*pp+vc.pts[v[k]]*0.5);
              vertices.push_back(*pp+vc.pts[v[k]+1]*0.5);
              counter++;
              k++;
            }
            else{
              printf(" ()");
            }
          } else {
            j=k+l;
            //printf("new Vector3(%gf,%gf,%gf),\r\n",*pp+vc.pts[v[k]-1]*0.5,pp[1]+vc.pts[v[k]]*0.5,pp[2]+vc.pts[v[k]+1]*0.5);
            vertices.push_back(*pp+vc.pts[v[k]-1]*0.5);
            vertices.push_back(*pp+vc.pts[v[k]]*0.5);
            vertices.push_back(*pp+vc.pts[v[k]+1]*0.5);
            counter++;
            k++;

            while(k<j){
              //printf("new Vector3(%gf,%gf,%gf),\r\n",*pp+vc.pts[v[k]-1]*0.5,pp[1]+vc.pts[v[k]]*0.5,pp[2]+vc.pts[v[k]+1]*0.5);
              vertices.push_back(*pp+vc.pts[v[k]-1]*0.5);
              vertices.push_back(*pp+vc.pts[v[k]]*0.5);
              vertices.push_back(*pp+vc.pts[v[k]+1]*0.5);
              counter++;
              k++;
            }

          }
        }
      }
    //}
      cells.push_back(counter++);*/
  }while(vl.inc());


  printf("done");
  con.draw_cells_gnuplot("random_points_v.gnu");
  con.print_custom("order=%o, vertices=%p","vertices.test");

  con.print_custom(
                "ID=%i, pos=(%x,%y,%z), vertices=%w, edges=%g, faces=%s, face_verts=%t"
                ,"face.test");
  printf("\n%lu,%lu\n", vertices.size(), cells.size());
  cellsPtr = new int[cells.size()];
  verticesPtr = new double[vertices.size()];
  std::copy(cells.begin(), cells.end(), (cellsPtr));
  std::copy(vertices.begin(), vertices.end(), (verticesPtr));
}




