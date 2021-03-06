// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/prevprime.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/prevprime.C
 * @brief NO DOC
 */
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <givaro/givintprime.h>
#include <givaro/givtimer.h>
#include <givaro/givinit.h>         // Givaro initialization


using namespace Givaro;




int main(int argc, char** argv)
{
//  Givaro::Init(&argc, &argv);


  IntPrimeDom IP;
  IntPrimeDom::Element m, ff;
  if (argc > 1) m = Integer(argv[1]);
  else std::cin >> m;
        Timer tim; tim.clear(); tim.start();
        IP.prevprimein(m);
        tim.stop();
        cout << m << endl;
        cerr << tim << endl;

//  Givaro::End();

  return 0;
}

