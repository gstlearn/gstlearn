/**
 * This test is meant to check the Environment
 */

#include "Basic/OptCst.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "Basic/File.hpp"

int main(int /*argc*/, char **/*argv*/)
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  // Management of OptDbg
  mestitle(0,"Checking 'OptDbg' functionalities");
  OptDbg::define(EDbg::NBGH);
  OptDbg::display();
  message("Current value of Neighborhood Debug = %d\n",
          OptDbg::query(EDbg::NBGH));

  // Management of OptCst
  mestitle(0,"Checking 'OptCst' functionalities");
  OptCst::display();

  VectorDouble vect = ut_vector_simulate_gaussian(20);
  message("\n");
  ut_vector_display("Vector Display (initial)", vect);

  OptCst::define(ECst::NTCAR,  12.);
  OptCst::define(ECst::NTDEC,   4.);
  OptCst::define(ECst::NTBATCH, 6.);
  OptCst::display();
  message("\n");
  ut_vector_display("Vector Display (modified)",vect);

  // Management of OptCustom
  mestitle(0,"Checking 'OptCustom' functionalities");
  OptCustom::define("My first Trial",  12.);
  OptCustom::define("My second Trial", 21.);
  OptCustom::define("My third Trial",  32.);
  message("Questioning known keyword = %lf\n",OptCustom::query("My first Trial"));
  message("Questioning unknown keyword (with default) = %lf\n",OptCustom::query("bidon",123.));
  message("Questioning unknown keyword (without default) = %lf\n",OptCustom::query("bidon"));
  OptCustom::undefine("Bidon"); // Does not do anything as keyword is not registered
  OptCustom::undefine("My second Trial");
  OptCustom::display();
}
