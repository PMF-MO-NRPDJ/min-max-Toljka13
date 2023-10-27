#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

#include "dune/common/parallel/mpihelper.hh"
#include <dune/common/exceptions.hh> 

#include <dune/grid/uggrid.hh>  
#include <dune/grid/common/gridinfo.hh> 
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

// Izračunaj kut između p1-p2 i p3-p2
template <typename Point>
double kut(Point const & p1, Point const & p2, Point const & p3)
{
   // računanje kuta (VAŠ KOD) ide ovdje    
   // Point će biti Dune::FieldVector<double, dim>. Norma se računa 
   // pomoću funkcije članice klase two_norm(), a skalarni produkt 
   // pomoću funkcije članice klase dot(). Pogledati dokumentaciju klase
   //  Dune::FieldVector<double, dim>.

    Dune::FieldVector<double, 2> v1 = p1 - p2; //vektor formiran od točaka p2 i p1
    Dune::FieldVector<double, 2> v2 = p3 - p2; //vektor formiran od točaka p2 i p3

    double skProd = dot(v1,v2); //skalarni produkt
    double prodNormi = v1.two_norm() * v2.two_norm(); //produkt normi

    double kut = std::acos(skProd / prodNormi); //cos(a,b) = (a*b) / |a||b|

    return 180*kut/M_PI;
}

int main(int argc, char** argv)
{
    const int dim = 2;
    using GridType = Dune::UGGrid<dim>;
    using LeafGridView = GridType::LeafGridView;
     
    // UČITATI 2D MREŽU IZ GMSH DATOTEKE KOJA JE ZADANA KAO ARGUMENT KOMANDNE LINIJE.
    
    if (argc == 1) {
            std::cerr << "Usage: " << argv[0] << " ime_grid_datoteke.msh"
                      << std::endl;
            std::exit(1);
        }
    Dune::MPIHelper::instance(argc, argv);

    bool verbosity = true;
    bool insertBoundarySegments = false;  // Bez toga Dune::GmshReader zna podbaciti (u 3D)

    std::unique_ptr<GridType> pgrid = Dune::GmshReader<GridType>::read(argv[1], verbosity, insertBoundarySegments);

    int no_r = std::stoi(argv[2]);  // broj profinjenja
    pgrid->globalRefine(no_r);     // profini mrežu

    auto gridView = pgrid->leafGridView();

    int brojac=0, i;
    double max = -1, min = 200, kut1;

    for(auto const & element : elements(gridView))
    {

/*     VAŠ KOD dolazi ovdje.
 *     RAČUNATI MIN I MAX KUT U SVAKOM ELEMENTU. 

*/

     for(i=0; i<=2; i++){
         kut1 = kut(element.geometry().corner(i%3),element.geometry().corner((i+1)%3),element.geometry().corner((i+2)%3));

         if(kut1 < min)
             min = kut1;
         if(kut1 > max)
             max = kut1;
     }

     brojac++;
    }

   // ISPISATI BROJ ELEMENATA; MINIMALNI I MAKSIMALNI KUT U STUPNJEVIMA:
    
    std::cout << "Broj elemenata mreže iznosi " << brojac << std::endl;
    std::cout << "Najveći kut u stupnjevima iznosi " << max << "°, a najmanji iznosi " << min << "°" << std::endl;

    // Ispis mreže u VTK formatu (u datoteku poluvijenac.vtu)
    Dune::VTKWriter<LeafGridView> vtkwriter(gridView);
    vtkwriter.write("poluvijenac");


    return 0;
}
