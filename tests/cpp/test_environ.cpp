/**
 * This test is meant to check the Environment
 */

#include "Basic/OptCst.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"

int main(int argc,char **argv)
{
  // Management of OptDbg
  mestitle(0,"Checking 'OptDbg' functionalities");
  OptDbg::define(EDbg::NBGH, true);
  OptDbg::display();
  message("Current value of Neighborhood Debug = %d\n",
          OptDbg::query(EDbg::NBGH));

  // Management of OptCst
  mestitle(0,"Checking 'OptCst' functionalities");
  VectorDouble vect = ut_vector_simulate_gaussian(20);

  OptCst::display();
  message("\n");
  ut_vector_display("Vector Display (initial)", vect);

  OptCst::define(ECst::NTCAR,  12.);
  OptCst::define(ECst::NTDEC,   4.);
  OptCst::define(ECst::NTBATCH, 6.);
  OptCst::display();
  message("\n");
  ut_vector_display("Vector Display (modified)",vect);
}
