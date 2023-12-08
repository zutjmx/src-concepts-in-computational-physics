// FileName.cpp : implementation file
//

#include "stdafx.h"
#include "FileName.h"
#include <direct.h>
#include <stdlib.h>

CString CFileName::m_sLoggedDir="";

// CFileName

/***************************************************************************/
int CFileName::SetLoggedDir()
/***************************************************************************/
{
   char* st;
   if ((st = _getcwd(NULL,0)) == NULL)          // get current directory
      return 0;
   m_sLoggedDir = st;
   m_sLoggedDir += '\\';
   free(st);
   return 1; 
}
/***************************************************************************/
CString CFileName::GetName()
/***************************************************************************/
{
    int pos = ReverseFind('\\');
    CString temp;
    if (pos == -1)
        temp.Empty();
    else {
        int len = GetLength()-pos-1;
        temp = Right(len);
    }
    return temp;
}
/***************************************************************************/
CString CFileName::GetNameSansExtension()
/***************************************************************************/
{
   int x = ReverseFind('\\'),
      pos = ReverseFind('.');
   CString temp;
   if (pos == -1 && x == -1)
      temp = *this;
   else {
      if (pos == -1)
         pos = GetLength();
      if (x == -1)
         x = 0;
      int len = pos-1-x;
      if (len)
         temp = Mid(x+1,len);
      else
         temp.Empty();       
    }
    return temp;
}
/***************************************************************************/
CString CFileName::GetExtension()
/***************************************************************************/
{
    int pos = ReverseFind('.');
    CString temp;
    if (pos == -1)
        temp.Empty();
    else {
        int len = GetLength()-pos-1;
        temp = Right(len);
    }
    return temp;
}
/***************************************************************************/
CString CFileName::GetPath()
/***************************************************************************/
{
    int pos = ReverseFind('\\');
    CString temp;
    if (pos == -1)
        temp.Empty();
    else
        temp = Left(pos+1);
    return temp;
}
/***************************************************************************/
void CFileName::ChangeFileExtension(const char *pszExt)
/***************************************************************************/
{
   int pos = ReverseFind('.');
   size_t x = GetLength()+4;
   LPTSTR buff = GetBuffer(static_cast<int>(x));
   if (pos == -1) {
      strcat_s(buff,x,".");
      strcat_s(buff,x,pszExt);
   } else {
      *(buff+pos+1) = '\0';
      strcat_s(buff,x,pszExt);
   }
   ReleaseBuffer();
}
/***************************************************************************/
void CFileName::ChangeFileExtension(UINT uExtID)
/***************************************************************************/
{
   int pos = ReverseFind('.');
   CString sExt;
   size_t x = GetLength()+4;
   sExt.LoadString(uExtID);
   LPTSTR buff = GetBuffer(static_cast<int>(x));
   if (pos == -1) {
      strcat_s(buff,x,".");
      strcat_s(buff,x,sExt);
   } else {
      *(buff+pos+1) = '\0';
      strcat_s(buff,x,sExt);
   }
   ReleaseBuffer();
}
