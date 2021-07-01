#ifndef APARAM_HPP
#define APARAM_HPP

class AParam
{
public:
  AParam(){}
  virtual ~AParam(){}

  virtual bool checkConsistence() const = 0;
private:
};
#endif
