prefix = input('prefix?', 's')
suffix = input('suffix?', 's')
idlist = input('ID list?', 's')
newfolder = input('new folder name?', 's')

idlist=load(idlist);
idlist=reshape(idlist,1,[]);
mkdir(newfolder);

for id_cnt=idlist
    idstr=num2str(id_cnt,'%5.5u')
    src_str = strcat(pwd,'/', prefix, idstr, suffix);
    dest_str = strcat(pwd,'/', newfolder, '/', prefix, idstr, suffix);
    copyfile(src_str, dest_str)
end