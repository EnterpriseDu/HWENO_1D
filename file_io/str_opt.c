#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>




/* This function examine whether a string
 * represents a real number.
 * Transform the string represents a
 * negtive number into a string represents
 * the opposite one and return its' sign.
 * It returns 0 if the string do not
 * represents a real number.
 * In the result given by this function, there will
 * be only one 'e' in the string, and the
 * only position for '-' is behind 'e', and
 * there can be only one dot in the string
 * and the only position for it is before 'e'.
 */
int format_string(char * str)
{
  int i = 0, length = 0, j = 0;
  int sign = 1;
  int flag_dot = 0; //the number of dots in the string
                    //should be at most one
  int pos_dot = 0;
  int flag_e = 0;
  int pos_e = 0;

  //printf("%send\n", str);

  length = strlen(str);

  for(j = 0; j < length; ++j)
  {
    if(!((str[j] == 46)||(str[j] == 45)||(str[j] == 43)||(str[j] == 69)||(str[j] == 101)||((str[j] > 47) && (str[j] < 58))))
      return 0;

    if((str[j] == 69) || (str[j] == 101))
    {
      str[j] = 101;
      flag_e += 1;
      pos_e = j;
    }
  }

  //there could not be more than one 'e' in one number
  if(flag_e > 1)
    return 0;
  if((flag_e) && (pos_e == 0))
    return 0;
  if((flag_e) && (pos_e == length-1))
    return 0;
  // a dot only could not be a number
  if((str[0] == 46) && (length == 1))
    return 0;
  // a '-' only could not be a number
  if(str[0] == 45)
  {
    if(length == 1)
      return 0;
    sign = -1;
  }
  // a '+' only could not be a number
  if(str[0] == 43)
    if(length == 1)
      return 0;

  // eliminate '-/+' from the string
  if((sign < 0) || (str[0] == 43))
  {
    for(i = 0; i < length; ++i)
      str[i] = str[i+1];
    length -= 1;
    pos_e -= 1;
    if(pos_e == 0)
      return 0;
  }

  if(flag_e)
  {
    for(i = 0; i < length; ++i)
    {
      if((str[i] == 45) || (str[i] == 43))
      {
      // after eliminating '-/+', the only possible position for '-/+'
      // is behind 'e'
        if((i-pos_e) != 1)
	  return 0;
	else if(i == length-1)
	  return 0;
      }
      //there could not be two dots in one number
      else if((str[i] == 46) && (flag_dot > 0))
        return 0;
      else if(str[i] == 46)
      {
        flag_dot += 1;
        pos_dot = i;
      }
    }
    if((flag_dot) && (pos_dot > (pos_e-2)))
      return 0;
  }
  else
  {
    for(i = 0; i < length; ++i)
    {
      // after eliminating '-/+', there should be no more '-/+'
      if((str[i] == 45) || (str[i] == 43))
        return 0;
      else if((str[i] == 46) && (flag_dot > 0))
        return 0;
      else if(str[i] == 46)
        flag_dot += 1;
    }
  }

  return sign;
}

/* this function transform a string
 * consisting '1', '2', ..., and '.'
 * into the real number it represent
 */
double str2num(char * number)
{
  double result = 0.0, super_script = 0.0;
  int idx = 0, dot = -2;
  int i = 0, j = 0, power, k = 0;
  int length = 0;
  int pos_e = 0;
  char * after_e = number;
  int sign = 1;

  length = strlen(number);

  for(j = 0; j < length; ++j)
    if(number[j] == 101)
      pos_e = j;

  /*
  for(k = 0; k < length; ++k)
    printf("%c", number[k]);
  printf("\t");
  */

  if(pos_e)
  {
    after_e = number + pos_e + 1;
    number[pos_e] = 0;
    result = str2num(number);
    if(after_e[0] == 45)
    {
      sign = -1;
      after_e += 1;
    }
    else if(after_e[0] == 43)
      after_e += 1;
    super_script = str2num(after_e);
    result = result * pow(10.0, sign * super_script);
  }
  else
  {
    while(number[idx] != 0)
    {
      if(number[idx] == 46)
      {
        dot = idx - 1;
	idx = 0;
	break;
      }
      ++idx;
    }

    if(dot == -2)
      dot = idx - 1;

    for (i = 0; i <= dot; ++i)
      result += (double)(number[i] - 48)*pow(10, dot - i);

    dot += 1;
    double addon = 0.0;
    for (i = 1; i < length - dot; ++i)
      result += (double)(number[dot+i] - 48)*pow(0.1, i);
  }
  
  return result;
}


/* this function counts how many lines are there
 * the initial data file.
 */
int data_pre_read_line(FILE * fp, int err_code)
{
  int line = 0,                 column = 0, M = 0;
  int flag = 0;  /* We need to know how many numbers are there in
		  * the initial data. "flag" helps us to count.
		  * We read characters one by one from the data
		  * file. The value of "flag" is 1 when read a
		  * number-using character (i.e. 0, 1, 2, and so
		  * on and the dot), while is 0 when read a 
		  * non-number-using character. 
		  */
  char ch;
  int num_in_line = 0;

  while((ch = getc(fp)) != EOF)
  {
    ++num_in_line;

    /* We get a empty charactor and the last charactor is a number-using one.
     * That means we just finished reading a string representing a number.
     */
    if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (flag))
    {
      ++column; //add the count of the column by 1
      flag = 0; //reset the value of flag to 0

      /* If the empty charactor wo get just now is a line-turning,
       * check if it has the same length as the previous ones does.
       */
      if(ch == '\n')
      {
	if(!line)
	  M = column;
	else if(column != M)
        {
	  printf("LINES NOT EQUAL, source=%d, line=%d, M=%d, column=%d\n", err_code, line, M, column);
	  return -err_code-1;
	}
	//if no problem, reset and go on
	++line;
	column = 0;
        num_in_line = 0;
      }
    }

    /* We get a number-using charactor and the last one is not.
     * That means we just begin reading a string representing a number.
     */
    else if( ((ch == 46)||(ch == 45)||(ch == 43)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (!flag) )
      flag = 1;

    /* We get an empty charactor and the last one is one, too.
     * nothing special, just check if it is a line turing.
     */
    else if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (!flag))
    {
      if(ch == '\n')
      {
	if(!line)
	  M = column;
	else if(column != M)
        {
	  printf("LINES NOT EQUAL, source=%d, line=%d, M=%d, column=%d\n", err_code, line, M, column);
	  return -err_code-1;
	}

	++line;
	column = 0;
        num_in_line = 0;
      }
    }

    /* We get a number-using charactor and the last one is one, too.
     * Nothing spectial.
     */
    else if( ((ch == 46)||(ch == 45)||(ch == 43)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (flag) )
      continue;

    /* The charactor we just get is neither an empty one nor a number-using one
     * Something wrong happened.
     */
    else
    {
      printf("Input (%d, %d, %d) contains illigal character(ASCII=%d, flag=%d)!\n", err_code, line+1, num_in_line, (int)ch, flag);
      return -err_code-9;
    }
  }



  if(flag)
    ++column;
  if(column)
    ++line;

  return line;
}


/*
 * This function tells us how many columns
 * are there in the file.
 */
int data_pre_read_column(FILE * fp, int err_code)
{
  int column = 0;
  int flag = 0;  /* We need to know how many numbers are there in
		  * the initial data. "flag" helps us to count.
		  * We read characters one by one from the data
		  * file. The value of "flag" is 1 when read a
		  * number-using character (i.e. 0, 1, 2, and so
		  * on and the dot), while is 0 when read a 
		  * non-number-using character. 
		  */
  char ch;
  int num_in_line = 0;

  while((ch = getc(fp)) != EOF)
  {
    ++num_in_line;
    if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (flag))
    {
      ++column;
      flag = 0;
      if(ch == '\n')
        break;
    }
    else if( ((ch == 46)||(ch == 45)||(ch == 43)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (!flag) )
      flag = 1;
    else if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (!flag))
    {
      if(ch == '\n')
	break;
      continue;
    }
    else if( ((ch == 46)||(ch == 45)||(ch == 43)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (flag) )
      continue;
    else
    {
      printf("Input (%d, %d) contains illigal character(ASCII=%d, flag=%d)!\n", err_code, num_in_line, (int)ch, flag);
      return -err_code-9;
    }
  }

  if(flag)
    ++column;

  return column;
}


/* This function reads the initial data file
 * to generate the initial data.
 * It returns 0 if successfully read the file,
 * while returns the index of the wrong entry.
 *
 * num = the expected number of entries in the file.
 */
int data_read(FILE * fp, double * U, int num)
{
  int idx = 0, j = 0; //j is a frequently used index for spatial variables
  char number[100];
  char ch;
  int sign = 1;

  while((ch = getc(fp)) != EOF)
  {
    if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (idx))
    {
      number[idx] = 0;
      idx = 0;

      sign = format_string(number);
      if(!sign)
	return j;
      // before EOF, j should always be smaller than num
      if(j == num)
	return j;

      U[j++] = sign * str2num(number);
    }
    else if((ch == 46)||(ch == 45)||(ch == 43)||((ch > 47) && (ch < 58))||(ch == 69)||(ch == 101))
      number[idx++] = ch;
  }

  // if the file is ended by a non-empty charactor before EOF, this block will be used
  if(idx)
  {
    sign = format_string(number);
    if(!sign)
      return j;
    if(j == num)
      return j;

    U[j] = sign * str2num(number);
  }

  return 0;
}



/* this function counts how many lines are there
 * the configuration file.
 */
int string_pre_read_line(FILE * fp, int err_code)
{
  int line = 0,                 column = 0, M = 0;
  int flag = 0;  /* We need to know how many numbers are there in
		  * the initial data. "flag" helps us to count.
		  * We read characters one by one from the data
		  * file. The value of "flag" is 1 when read a
		  * number-using character (i.e. 0, 1, 2, and so
		  * on and the dot), while is 0 when read a 
		  * non-number-using character. 
		  */
  char ch;
  int num_in_line = 0;

  while((ch = getc(fp)) != EOF)
  {
    ++num_in_line;

    /* We get a empty charactor and the last charactor is a number-using one.
     * That means we just finished reading a string representing a number.
     */
    if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (flag))
    {
      ++column; //add the count of the column by 1
      flag = 0; //reset the value of flag to 0

      /* If the empty charactor wo get just now is a line-turning,
       * check if it has the same length as the previous ones does.
       */
      if(ch == '\n')
      {
	if(!line)
	  M = column;
	else if(column != M)
        {
	  printf("LINES NOT EQUAL, source=%d, line=%d, M=%d, column=%d\n", err_code, line, M, column);
	  return -err_code-1;
	}
	//if no problem, reset and go on
	++line;
	column = 0;
        num_in_line = 0;
      }
    }

    /* We get a number-using charactor and the last one is not.
     * That means we just begin reading a string representing a number.
     */
    else if( ((ch == 46)||(ch == 45)||(ch == 43)||((ch > 64)&&(ch < 91))||((ch > 96)&&(ch < 123))||((ch > 47)&&(ch < 58))) && (!flag) )
      flag = 1;

    /* We get an empty charactor and the last one is one, too.
     * nothing special, just check if it is a line turing.
     */
    else if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (!flag))
    {
      if(ch == '\n')
      {
	if(!line)
	  M = column;
	else if(column != M)
        {
	  printf("LINES NOT EQUAL, source=%d, line=%d, M=%d, column=%d\n", err_code, line, M, column);
	  return -err_code-1;
	}

	++line;
	column = 0;
        num_in_line = 0;
      }
    }

    /* We get a number-using charactor and the last one is one, too.
     * Nothing spectial.
     */
    else if( ((ch == 46)||(ch == 45)||(ch == 43)||((ch > 64)&&(ch < 91))||((ch > 96)&&(ch < 123))||((ch > 47)&&(ch < 58))) && (flag) )
      continue;

    /* The charactor we just get is neither an empty one nor a number-using one
     * Something wrong happened.
     */
    else
    {
      printf("Input (%d, %d, %d) contains illigal character(ASCII=%d, flag=%d)!\n", err_code, line+1, num_in_line, (int)ch, flag);
      return -err_code-9;
    }
  }



  if(flag)
    ++column;
  if(column)
  {
    if(column != M)
    {
      if(ch == EOF)
	printf("EOF\n");
      printf("LINES NOT EQUAL, source=%d, line=%d, M=%d, column=%d\n", err_code, line, M, column);
      return -err_code-1;
    }
    ++line;
  }

  return line;
}


/*
 * This function tells us how many columns
 * are there in the file.
 */
int string_pre_read_column(FILE * fp, int err_code)
{
  int column = 0;
  int flag = 0;  /* We need to know how many numbers are there in
		  * the initial data. "flag" helps us to count.
		  * We read characters one by one from the data
		  * file. The value of "flag" is 1 when read a
		  * number-using character (i.e. 0, 1, 2, and so
		  * on and the dot), while is 0 when read a 
		  * non-number-using character. 
		  */
  char ch;
  int num_in_line = 0;

  while((ch = getc(fp)) != EOF)
  {
    ++num_in_line;
    if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (flag))
    {
      ++column;
      flag = 0;
      if(ch == '\n')
        break;
    }
    else if( ((ch == 46)||(ch == 45)||(ch == 43)||((ch > 64)&&(ch < 91))||((ch > 96)&&(ch < 123))||((ch > 47)&&(ch < 58))) && (!flag) )
      flag = 1;
    else if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (!flag))
    {
      if(ch == '\n')
	break;
      continue;
    }
    else if( ((ch == 46)||(ch == 45)||(ch == 43)||((ch > 64)&&(ch < 91))||((ch > 96)&&(ch < 123))||((ch > 47)&&(ch < 58))) && (flag) )
      continue;
    else
    {
      printf("Input (%d, %d) contains illigal character(ASCII=%d, flag=%d)!\n", err_code, num_in_line, (int)ch, flag);
      return -err_code-9;
    }
  }

  if(flag)
    ++column;

  return column;
}



int string_read(FILE * fp, char * str_conf[], int const n_conf, int const L)
{
  int idx = 0, k = 0;
  char ch;

  while((ch = getc(fp)) != EOF)
  {
    if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (idx))
    {
      str_conf[k][idx] = 0;
      ++k;
      idx = 0;
    }
    else if((ch == 46)||(ch == 45)||(ch == 43)||((ch > 64)&&(ch < 91))||((ch > 96)&&(ch < 123))||((ch > 47)&&(ch < 58)))
    {
      if(k == n_conf)
	return k;
      if(idx == L)
	return -k;

      str_conf[k][idx++] = ch;
    }
  }


  return 0;
}
