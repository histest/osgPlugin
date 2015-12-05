
#include "stdafx.h"
#include <stdio.h>
#include <string.h>

int	ReadFileLine(FILE *inFile, char chrLineBuffer[82])
{
	char chrIgnore[40];
	char *chrCode;
	int	i;
	int	iStatus = 0;

	chrCode = fgets(chrLineBuffer,82,inFile);

	/*一行文本超过81个字符，文件指针转到下一行*/
	if ( (strlen(chrLineBuffer) == 81) && (chrLineBuffer[80] != '\n') )
	{
		fgets(chrIgnore,40,inFile);                
		chrLineBuffer[81] = '\0';
	}

	/*转换Fortran的双精度标志'D'为C的'E'*/
	for (i=0;i<82;i++)
	{
		if (chrLineBuffer[i] == '\0') break;
		if (chrLineBuffer[i] ==  'D') chrLineBuffer[i] = 'E';
	}

	if ( chrCode == NULL )	iStatus = EOF;
	return iStatus;
}
