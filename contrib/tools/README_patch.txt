The patch.exe in this directory was created using the original patch.exe binary from 
http://gnuwin32.sourceforge.net/packages/patch.htm with its code remaining unmodified.
The only difference is that it contains an embedded manifest, allowing this patch.exe to run without admin privileges.



Technical details:

Windows in its infinite wisdom (based on a filename heuristic) thinks that the original GNU patch.exe requires elevation (i.e. admin privileges)
-- which of course is not the case -- but Windows still insists.
Thus, you either have to 

A) Run it as an admin (if your account allows and bend to the will of Windows)

or

B) Supply an XML-based manifest file (see below) named patch.exe.manifest with resides
side-by-side with the patch.exe AND is older(!) in terms of modification date.
For this, you might need to run touch.exe (e.g. http://www.binarez.com/touch_dot_exe/#download) on patch.exe.

or

C) Embed the manifest into patch.exe using MT.exe (which comes with older versions of Visual Studio or Windows SDK) by running
   
  > mt -manifest patch.exe.manifest -outputresource:patch.exe;1

The modified patch.exe with an embedded manifest from solution C) is provided here.


---- CONTENT OF 'patch.exe.manifest'  -----------

<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<assembly xmlns="urn:schemas-microsoft-com:asm.v1" manifestVersion="1.0">
  <trustInfo xmlns="urn:schemas-microsoft-com:asm.v3">
    <security>
      <requestedPrivileges>
        <requestedExecutionLevel level="asInvoker" uiAccess="false"/>
      </requestedPrivileges>
    </security>
  </trustInfo>
</assembly>