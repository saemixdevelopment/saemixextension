function addmenumlx(hObject)
%ADDMENUMLX  add menu to figure, called by MONOLIX

%  Marc Lavielle
%  Version 1.1 ;  2004/02/18

nh=length(hObject);
for k=1:nh
set(hObject(k),'menu','none')
uimenu(hObject(k),'Label','Save','Callback','[f,p]=uiputfile(''*.fig'',''Save as''); saveas(gcf,[p f]);');
uimenu(hObject(k),'Label','Print','Callback','printdlg(''-crossplatform'')');
uimenu(hObject(k),'Label','Close','Callback','close');
uimenu(hObject(k),'Label','Zoom on','Callback','zoom on');
uimenu(hObject(k),'Label','Zoom out','Callback','zoom out');
end;

