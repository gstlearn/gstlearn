/******************************************************************************/
/*                                                                            */
/*                            VB C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: VB Team                                                     */
/* Website: https://VB.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/ASerializable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "DataBase/DbData.hpp"
#include "DataBase/Dictionary.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of DbData
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_DataBase-");

  DbData data{};
  data.printContents();

  // Checking the different types of Columns that can be added to a DbData
  mestitle(1, "Adding Columns to a DbData");
  auto nech = 3;
  data.addColumn("VD", VectorDouble{1., 2., 3.}, RoleID{ERole::X});
  data.addColumn("VI", VectorInt{5, 6, 7}, RoleID{ERole::Z});
  data.addColumn("VS", VectorString{"foo", "bar", "baz"}, RoleID{ERole::Z, 1});
  data.addColumn("VB", VectorBool{true, false, true});
  data.addColumn("VDS", VH::sequenceVD(2., 19., 1.), RoleID(ERole::F), 6);
  data.addColumn("VIS", VH::sequence(15, 3, 2), RoleID(ERole::Z, 2), 5);
  data.printContents();

  // Checking the addition of a Categorical variable
  message("Creating the Dictionary for the categorical variable\n");
  VectorInt keys{1, 3, 7};
  VectorString labels{"red", "blue", "green"};
  Dictionary dict;
  for (size_t i = 0; i < keys.size(); ++i) dict.addCategory(keys[i], labels[i]);
  auto ncategory = dict.getNCategories();

  message("Creating the VectorCategory for the categorical variable\n");
  VectorCategory vc(nech, dict);
  for (Id i = 0; i < nech; ++i)
  {
    auto j = static_cast<Id>(i % ncategory);
    vc.setCategory(i, keys[j]);
  }
  data.addColumn("VC", VectorCategory(vc));
  data.printContents();

  // Checking the use of Neutral Files
  mestitle(1, "Saving and recovering a DbData from a Neutral File");
  data.printContents("Before Saving in a Neutral File");
  std::cout << "VD: " << data.getColumn<VectorDouble>("VD");
  std::cout << "VI: " << data.getColumn<VectorInt>("VI");
  std::cout << "VS: " << data.getColumn<VectorString>("VS");
  std::cout << "VB: " << data.getColumn<VectorBool>("VB");
  std::cout << "VDS: " << data.getColumn<VectorDouble>("VDS");
  std::cout << "VIS: " << data.getColumn<VectorInt>("VIS");
  std::cout << "VC: " << data.getColumn<VectorCategory>("VC");
  data.dumpToNF("test_DataBase.NF");

  message("\nAfter recovering from the Neutral File\n");
  auto* data2 = DbData::createFromNF("test_DataBase.NF", true);
  std::cout << "VD: " << data2->getColumn<VectorDouble>("VD");
  std::cout << "VI: " << data2->getColumn<VectorInt>("VI");
  std::cout << "VS: " << data2->getColumn<VectorString>("VS");
  std::cout << "VB: " << data2->getColumn<VectorBool>("VB");
  std::cout << "VDS: " << data2->getColumn<VectorDouble>("VDS");
  std::cout << "VIS: " << data2->getColumn<VectorInt>("VIS");
  std::cout << "VC: " << data2->getColumn<VectorCategory>("VC");
  delete data2;

  // Particular test for the case of Categorical variable
  mestitle(1, "Checking the use of a Categorical Variable");
  // Reading the category of the first sample and store it in 'cat'
  Category cat = data.getValue<Category>("VC", 0).value_or(getNA<Category>());
  // Setting the category 'cat' to the second sample
  data.setValue("VC", 1, cat);
  // Cheking the result on the dump of the whole vector
  std::cout << "VC: " << data.getColumn<VectorCategory>("VC");

  // Checking aliases for accessing values in a DbData
  mestitle(1, "Checking aliases");
  message("Various ways to get value of 'VIS' for version 0 at sample 2\n");
  auto is = 2;
  auto ind = 2;
  auto iv = 0;
  message("- By Name: %d\n", data.getValue<Id>("VIS", is).value_or(-1));
  message(
    "- By Name and Version: %d\n",
    data.getValue<Id>({"VIS", iv}, is).value_or(-1));
  message(
    "- Role: %d\n", data.getValue<Id>(RoleID{ERole::Z, ind}, is).value_or(-1));

  // Checking the different manners to refer to a Column in a DbData
  mestitle(1, "Different manners to refer to a Column");
  message("- by Name: %s\n", data.getName("VD").c_str());
  message("- by Name and Version: %s\n", data.getName({"VDS", 1}).c_str());
  message("- by Index: %s\n", data.getName(4).c_str());
  message("- by Index and Version: %s\n", data.getName({4, 1}).c_str());
  message("- by RoleID: %s\n", data.getName(RoleID(ERole::F)).c_str());
  message(
    "- by RoleID and Version: %s\n",
    data.getName({RoleID(ERole::F), 0}).c_str());
  message("- by Role: %s\n", data.getName(ERole::Z).c_str());

  // Retreiving values from a DbData and modifying them
  mestitle(1, "Retrieving and modifying values in a DbData");
  auto isample = 2;

  message("Initial values\n");
  message("%lf\n", data.getValue<double>("VD", isample).value_or(-1));
  message("%d\n", data.getValue<int>("VI", isample).value_or(-1));
  message(
    "%s\n", data.getValue<String>("VS", isample).value_or(String()).c_str());

  data.setValue("VD", isample, 4.);
  data.setValue("VI", isample, static_cast<Id>(8));
  data.setValue("VS", isample, String("mystring"));

  // Checking implicit conversions
  message("\nValues after modification\n");
  message("%lf\n", data.getValue<double>("VD", isample).value_or(-1));
  message("%d\n", data.getValue<int>("VI", isample).value_or(-1));
  message(
    "%s\n", data.getValue<String>("VS", isample).value_or(String()).c_str());

  mestitle(1, "Testing implicit conversions");

  double vd = data.getValue<double>("VD", 1).value_or(-1.);
  int vi = data.getValue<int>("VD", 1).value_or(-1);

  message("double from double: %lf\n", vd);
  message("int from double: %d\n", vi);

  double bad = data.getValue<double>("VS", 1).value_or(-999.);
  message("double from String: %lf\n", bad);

  // Checking proxy accesses
  mestitle(1, "Testing proxy access syntax");

  data.X()[2] = 123.;
  message(
    "Value of X(0) at sample 2: %lf\n",
    data.getValue<double>(RoleID{ERole::X}, 2).value_or(-1));

  data.F()(2, 3) = 12.34;
  message(
    "Value of F(0) at sample index 2 version index 3: %lf\n",
    data.getValue<double>(ColID(RoleID(ERole::F), 3), 2).value_or(-1));

  data.Z(2)(2, 4) = 999;
  message(
    "Value of Z(2) at sample index 2 version index 4: %d\n",
    data.getValue<int>(ColID(RoleID(ERole::Z, 2), 4), 2).value_or(-1));

  // Checking the aliases using ALL (or _) for accessing vectors of values in a DbData
  message("\nChecking the aliases using ALL(_)\n");
  data.printContents();
  auto FAll0 = data.F()();
  FAll0.dump("data.F()()");
  auto FAll1 = data.F()(_, _);
  FAll1.dump("data.F()(_, _)");
  auto FAll2 = data.F()(_);
  FAll2.dump("data.F()(_)");
  VectorDouble FAll3 = data.F()(_, 1);
  FAll3.dump("data.F()(_, 1)");
  auto FAll4 = data.F()(1, _);
  FAll4.dump("data.F()(1, _)");

  data.F() = 1234.;
  data.F()().dump("After global modification of F()");
  data.F()(2, _) = 5678.;
  data.F()().dump("After modification of F()(2, _)");
  data.F()(_, 3) = 91011.;
  data.F()().dump("After modification of F()(_, 3)");

  // Checking Columns with same Name and/or Role
  mestitle(1, "Adding Columns with the same Name and/or Role");
  data.printContents("Initial");
  data.addColumn("VD2", VectorDouble{1., 2., 3.}, RoleID{ERole::X, 10});
  data.addColumn("VI2", VectorDouble{1., 2., 3.}, RoleID{ERole::X, 0});
  data.addColumn("VI2", VectorDouble{1., 2., 3.}, RoleID{ERole::X, 10});
  data.printContents("\nFinal");

  // Provoking errors
  mestitle(1, "Modifying a Column with the wrong type");
  data.getColumn<VectorDouble>("VD").dump("Initial Column");
  data.setValue("VD", isample, 11123.);
  data.getColumn<VectorDouble>("VD").dump("After valid change");
  message("Trying to modify a Column with the wrong type\n");
  data.setValue("VD", isample, "String");

  mestitle(1, "Retrieving series of columns for various criteria");
  data.deleteAllColumns();
  data.addColumn("VD", VectorDouble{1., 2., 3.}, RoleID{ERole::X});
  data.addColumn("VDbis", VectorInt{5, 6, 7}, RoleID{ERole::X, 1});
  data.addColumn("VDback", VectorDouble{5, 6, 7}, RoleID{ERole::X, 2});
  data.printContents();

  message("\nMatching the name 'VDb*', the following columns are found:\n");
  auto colIDs = data.getColIDs("VDb*");
  for (auto& colID: colIDs)
    message("Matching Column: %s\n", colID.getDescr().c_str());

  message("\nMatching the role 'X', the following columns are found:\n");
  colIDs = data.getColIDs(ERole::X);
  for (auto& colID: colIDs)
    message("Matching Column: %s\n", colID.getDescr().c_str());

  return 0;
}
