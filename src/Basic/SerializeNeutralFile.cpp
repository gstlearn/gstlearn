/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/

#include "Basic/ASerializable.hpp"
#include "Basic/SerializeNeutralFile.hpp"

bool SerializeNeutralFile::fileOpenWrite(const ASerializable& parent,
                                         const String& filename,
                                         std::ofstream& os,
                                         bool verbose)
{
  // Close the stream if opened
  if (os.is_open()) os.close();
  // Build the multi-platform filename
  String filepath = ASerializable::buildFileName(2, filename, true);
  // Open new stream
  os.open(filepath, std::ios::out | std::ios::trunc);
  if (!os.is_open())
  {
    if (verbose) messerr("Error while opening %s", filepath.c_str());
    return false;
  }
  // Write the file type (class name)
  os << parent._getNFName() << std::endl;
  return os.good();
}

bool SerializeNeutralFile::fileOpenRead(const ASerializable& parent,
                                        const String& filename,
                                        std::ifstream& is,
                                        bool verbose)
{
  // Close the stream if opened
  if (is.is_open()) is.close();
  // Build the multi-platform filename
  String filepath = ASerializable::buildFileName(1, filename, true);
  // Open new stream
  is.open(filepath, std::ios::in);
  if (!is.is_open())
  {
    if (verbose) messerr("Error while opening %s", filepath.c_str());
    return false;
  }
  // Read and check the file type (class name)
  String type;
  is >> type;
  if (type != parent._getNFName())
  {
    if (verbose)
      messerr("The file %s has the wrong type (read: %s, expected: %s)", filepath.c_str(),
              type.c_str(), parent._getNFName().c_str());
    is.close();
    return false;
  }
  return is.good(); // Cannot be "end of file" already
}

bool SerializeNeutralFile::commentWrite(std::ostream& os, const String& comment)
{
  if (os.good())
  {
    if (comment.empty())
      os << std::endl;
    else
      os << "# " << comment << std::endl;
  }
  return os.good();
}

bool SerializeNeutralFile::tableWrite(std::ostream& os,
                                      const String& string,
                                      int ntab,
                                      const VectorDouble& tab)
{
  bool ret = true;
  VectorDouble loctab(ntab);
  for (int i = 0; i < ntab; i++) loctab[i] = tab[i];
  ret = ret && recordWriteVec<double>(os, string, loctab);
  return ret;
}

bool SerializeNeutralFile::tableRead(std::istream& is,
                                     const String& string,
                                     int ntab,
                                     double* tab)
{
  bool ret = true;
  VectorDouble loctab(ntab);
  ret = ret && recordReadVec<double>(is, string, loctab, ntab);
  if (!ret) return 1;
  for (int i = 0; i < ntab; i++) tab[i] = loctab[i];
  return ret;
}

bool SerializeNeutralFile::onlyBlanks(char* string)
{
  int number = static_cast<int>(strlen(string));
  for (int i = 0; i < number; i++)
  {
    if (string[i] != ' ') return false;
  }
  return true;
}
