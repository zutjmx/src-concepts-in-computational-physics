#pragma once


// CFileName

class CFileName : public CString
{

public:
   CFileName(void) {}
   CFileName(const char *str): CString(str) {}
   CFileName(CString str): CString(str) {}
   CString GetName(void);
   CString GetNameSansExtension(void);
   CString GetExtension(void);
   CString GetPath(void);
   void ChangeFileExtension(const char *);
   void ChangeFileExtension(UINT);
   const CFileName& operator=(const char *str)
      { CString::operator=(str); return *this; }
   const CFileName& operator=(CString str)
      { CString::operator=(str); return *this; }
   static int SetLoggedDir(void);
   static CString &GetLoggedDir(void) { return m_sLoggedDir; }

private:
   static CString m_sLoggedDir;        //  logged directory 
};


