
/*                       
	 Copyright (C) 2005 Tom Drummond

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public
     License as published by the Free Software Foundation; either
     version 2.1 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public
     License along with this library; if not, write to the Free Software
     Foundation, Inc.
     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/
#ifndef __MEMBASE_H
#define __MEMBASE_H


// the class *is* the data (so a local variable places it on the stack)
template<int Size>
class Stack {
protected:
  double my_values[Size];
#ifndef WIN32
} __attribute__ ((aligned(16)));
#else
};
#endif

// the class allocates and deallocates the data on the heap
template<int Size>
class Heap {
public:
  inline Heap(){my_values = new double[Size];}
  inline Heap(const Heap& copyof){
    my_values = new double[Size];
    memcpy(my_values,copyof.my_values,Size*sizeof(double));
  }

  inline Heap& operator=(const Heap& copyof){
    memcpy(my_values,copyof.my_values,Size*sizeof(double));
  }

  inline ~Heap(){delete [] my_values;}
protected:
  double* my_values;
};







#endif
