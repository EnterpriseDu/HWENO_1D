#include <stdlib.h>
#include <stdio.h>

#include "file_io.h"



int main()
{
  runList runhist;

  init_runList(&runhist);

  int j;

  
  for(j = 0; j < 4; ++j)
  {
    insert_runList(&runhist);
    locate_runList(j+1, &runhist);
    runhist.current->RcstrState[0] = 48+j;
  }

  runhist.current = runhist.head;
  for(j = 0; j < runhist.length; ++j)
  {
    printf("%c\n", runhist.current->RcstrState[0]);
    runhist.current = runhist.current->next;
  }
  printf("%p\n\n", runhist.current);

  delete_runList(&runhist);
  printf("%p\n", runhist.head);
  printf("%p\n", runhist.current);
  printf("%p\n\n", runhist.tail);

  printf("%d\n\n", check_runList(&runhist));
}
