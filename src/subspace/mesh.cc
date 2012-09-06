#include <iostream>
#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>

#include "subspace/mesh.hh"


#ifdef WINDOWS
  #include <direct.h> 
  #define GetCurrentDir _getcwd
  #define PathSeparator 92 //backslash
#else
  #include <unistd.h>
  #define GetCurrentDir getcwd
  #define PathSeparator 47 //forward slash
#endif


namespace subspace {

  void TriMesh::read(std::string file, std::string type){

    char cCurrentPath[FILENAME_MAX];
    if (!GetCurrentDir(cCurrentPath, FILENAME_MAX))
      {
	exit(1);
      }

    std::cout << "Current Dir:" << cCurrentPath << std::endl;
    
    //file.insert(file.begin(),(char) PathSeparator);
    //file.insert(0,cCurrentPath);

    if (type.compare("off")==0){
      std::ifstream mesh_Fin;
      file.append("."); file.append("off");    
      mesh_Fin.open(file.c_str());
      if (!mesh_Fin.is_open()){
	std::cerr << "Error: Cannot open file " << file << std::endl;
	exit(1);
      }

      mesh_Fin >> P;// mesh read
      if (!mesh_Fin){
	std::cerr << file.c_str() << " is not a manifold mesh" <<std::endl;
	exit(1); 
      }

      mesh_Fin.close();
    }


    if (!P.is_pure_triangle()) {std::cerr<< "Error: input mesh has non-triangle faces!" <<std::endl; exit(1);}
    if (!P.is_closed()) {std::cout<< "Warning: input mesh seems not to be watertight" <<std::endl;}

    int rm_component = P.keep_largest_connected_components(1);
    if (rm_component != 0) {std::cout<<"Warning: multiple shells detected, remove "<<rm_component<<" shells to keep the largest one"<<std::endl;}
  };


  void TriMesh::write(std::string file, std::string type){

    if (type.compare("off")==0){
      std::ofstream mesh_Fout;
      file.append("."); file.append(type);
      mesh_Fout.open(file.c_str());
      mesh_Fout << P;// mesh output
      std::cout << "Export mesh to: " << file <<  std::endl;
      mesh_Fout.close();
    }
  };

  void TriMesh::init_index(){
    //int n=(P.size_of_halfedges()+P.size_of_border_edges())/2;

    halfedge_num = P.size_of_halfedges();
    vertex_num = P.size_of_vertices();
    facet_num = P.size_of_facets();

    

    IH =  ISHalfedgeList(halfedge_num);
    IV =  ISVertexList(vertex_num);
    IF =  ISFacetList(facet_num);

    halfedge_vec.resize(halfedge_num);
    vertex_norm.resize(vertex_num);
    facet_norm .resize(facet_num);
    vertex_area.resize(vertex_num);
    facet_area.resize(facet_num);

    vertex_avg_len.resize(vertex_num);


    int index_count=0;
    for(Vertex_iterator vitr= P.vertices_begin();vitr!= P.vertices_end();
	IV[index_count]=vitr, vitr->index = index_count++, vitr++);
    index_count=0;
    for(Halfedge_iterator eitr= P.halfedges_begin();eitr!= P.halfedges_end();
	IH[index_count]=eitr, eitr->index = index_count++, eitr++);
    index_count=0;
    for(Facet_iterator fitr= P.facets_begin(); fitr!= P.facets_end(); 
	IF[index_count]=fitr, fitr->index = index_count++, fitr++);
  };


  void TriMesh::update_internal(){
    init_index();

    for (int i=0;i<halfedge_num;i++) 
      {
	Halfedge_handle h = IH[i];
	halfedge_vec[i] = h->vertex()->point() - h->prev()->vertex()->point();
      }

    for (int i=0;i<facet_num;i++)
      {
	Halfedge_handle h=IF[i]->halfedge();
	Vector normal = CGAL::cross_product(halfedge_vec[h->index], halfedge_vec[h->next()->index]);
	facet_area[i] = CGAL::sqrt(normal * normal);
	facet_norm[i] = normal / facet_area[i]; 
	facet_area[i] /= 2.;
      }

    for (int i=0; i<vertex_num; i++){
      HV_circulator hv=IV[i]->vertex_begin();
      double area = 0, total_len =0;
      Vector tmp;
      int k=0;

      do {
	if (hv->facet()==NULL) continue;
	area += facet_area[hv->facet()->index];

	tmp = halfedge_vec[hv->index];
	total_len += CGAL::sqrt(tmp * tmp);		    
	++k;      

      }while (++hv!=IV[i]->vertex_begin());
      vertex_area[i] = area / k;
      vertex_avg_len[i] = total_len / k;

    }

    facet2vertex_point_average( facet_norm, vertex_norm, Vector(0,0,0));
    for (int i=0;i<vertex_num;i++){
      vertex_norm[i] = vertex_norm[i] / CGAL::sqrt(vertex_norm[i] * vertex_norm[i]);
    }

    vertex_array= new double [3 * vertex_num];
    normal_array= new double [3 * vertex_num];
    facet_array = new unsigned int[3 * facet_num];

    for (int i=0; i<vertex_num; i++) {
      Point p=IV[i]->point();
      Vector n=vertex_norm[i];
      vertex_array[3*i] = p.x();
      vertex_array[3*i+1] = p.y();
      vertex_array[3*i+2] = p.z();
      normal_array[3*i] = n.x();
      normal_array[3*i+1] = n.y();
      normal_array[3*i+2] = n.z();
    }

    for (int i=0; i<facet_num; i++) {
      HF_circulator hc = IF[i]->facet_begin();
      int j=0;//preconditon: triangular mesh
      do{
	facet_array[3*i+j] = hc->vertex()->index;
	++j;
      }while(++hc != IF[i]->facet_begin());        
    }


  }

  void TriMesh::destroy(){
    delete [] vertex_array;
    delete [] normal_array;
    delete [] facet_array;
  }


}
