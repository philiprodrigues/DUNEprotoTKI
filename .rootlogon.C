{
  //need to put the follwoing line in ~/.rootrc
  //Rint.Logon: /some_path/rootlogon.C
  //or
  //https://root-forum.cern.ch/t/how-to-change-default-rootlogon-c-file/23883
  /*
    The order in which the rootlogon.C is look for is

        $PWD/.rootlogon.C
        $HOME/.rootlogon.C
        $ROOTETCDIR/system.rootlogon.C
  */
  const TString pwd = gSystem->pwd();
  cout<<endl;
  cout<<"Welcome to protoDUNETKI repository!"<<endl;
  cout<<"In path: "<<pwd<<endl;
  cout<<endl;

  const TString nugas=gSystem->Getenv("DUNEPROTOTKI");
  if(nugas!=""){
    const TString libs[]={"EG", "style"};
    for(unsigned int ii=0; ii<sizeof(libs)/sizeof(TString); ii++){
      printf("Loading %s ... %d\n",libs[ii].Data(), gSystem->Load("lib"+libs[ii]+".so"));
    }
  }
}
