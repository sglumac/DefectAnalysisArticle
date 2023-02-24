#include "article.h"

#include <string>

int main(int argc, char *argv[])
{
  if (argc <= 1)
  {
    fmi::examples::article::results_for_all_figures();
  }
  else
  {
    int figure = std::stoi(argv[1]);
    switch(figure)
    {
      case 1:
        fmi::examples::article::experiment_1();
        break;
      case 2:
        fmi::examples::article::experiment_2();
        break;
      case 3:
        fmi::examples::article::experiment_3();
        break;
      case 4:
        fmi::examples::article::experiment_4();
        break;
      case 5:
        fmi::examples::article::experiment_5();
        break;
      case 6:
        fmi::examples::article::experiment_6();
        break;
    }
  }
  
  return 0;
}
