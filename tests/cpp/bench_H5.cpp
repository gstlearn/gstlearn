#include <Basic/Timer.hpp>
#include <Db/DbGrid.hpp>

int test_NF(const DbGrid& db)
{
  Timer tm;
  db.dumpToNF("dbgrid.nf");
  auto db2 = std::unique_ptr<DbGrid>(DbGrid::createFromNF("dbgrid.nf"));
  if (db2 != nullptr)
  {
    db2->dumpToNF("dbgrid2.nf");
  }
  else
  {
    messerr("Cannot deserialize `dbgrid.nf'");
    return 1;
  }
  tm.displayIntervalMilliseconds("Serialize + Deserialize + Serialize Neutral File",
                                 2500);
  return 0;
}

int test_HDF5(const DbGrid& db)
{
  Timer tm;
  db.dumpToH5("dbgrid.h5");
  auto db2 = std::unique_ptr<DbGrid>(DbGrid::createFromH5("dbgrid.h5"));
  if (db2 != nullptr)
  {
    db2->dumpToH5("dbgrid2.h5");
  }
  else
  {
    messerr("Cannot deserialize `dbgrid.h5'");
    return 1;
  }
  tm.displayIntervalMilliseconds("Serialize + Deserialize + Serialize HDF5", 80);

  return 0;
}

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  const VectorInt dims {100, 100, 100};
  auto db = std::unique_ptr<DbGrid>(DbGrid::create(dims));
  int ret {};
  ret += test_NF(*db);
  ret += test_HDF5(*db);
  return ret;
}
