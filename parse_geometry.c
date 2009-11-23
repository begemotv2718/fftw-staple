int parse_staple( char **rest, char *result){
  char *geometry;
  geometry=*rest;
  if(geometry[0]=='t' || geometry[0]=='b' || geometry[0]=='l' || geometry[0]=='r'){
    *result=geometry[0];
    (*rest)++;
    return(1);
  }
  return(0);
}
int parse_int(char **rest, int *result){
  *result=0;
  if(!((**rest>='0')&&(**rest<='9')))
    return(0);
  while((**rest>='0')&&(**rest<='9')){
    *result*=10;
    *result+=**rest-'0';
    (*rest)++;
  }
  return(1);
}
int parse_x( char **rest){
  if(*rest[0]=='x'){
    (*rest)++;
    return(1);
  }
  return(0);
}  
int parse_geometry(char *geometry, int *size_x, int *size_y, char *staple){
  char *rest;
  rest=geometry;
  return(
  parse_staple(&rest,staple) &&
  parse_int(&rest,size_x) &&
  parse_x(&rest) &&
  parse_int(&rest,size_y));
}

