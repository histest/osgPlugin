
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

	/*һ���ı�����81���ַ����ļ�ָ��ת����һ��*/
	if ( (strlen(chrLineBuffer) == 81) && (chrLineBuffer[80] != '\n') )
	{
		fgets(chrIgnore,40,inFile);                
		chrLineBuffer[81] = '\0';
	}

	/*ת��Fortran��˫���ȱ�־'D'ΪC��'E'*/
	for (i=0;i<82;i++)
	{
		if (chrLineBuffer[i] == '\0') break;
		if (chrLineBuffer[i] ==  'D') chrLineBuffer[i] = 'E';
	}

	if ( chrCode == NULL )	iStatus = EOF;
	return iStatus;
}
