#include <stdio.h>
#include <stdlib.h>
#include <Imlib2.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <math.h>

#define SIZE_X 100
#define SIZE_Y 500
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))


typedef struct{
  int x1,y1,x2,y2;
} rectangle;

void make_rect(int x1, int y1, int x2, int y2, rectangle *rect){
  rect->x1=x1;
  rect->y1=y1;
  rect->x2=x2;
  rect->y2=y2;
}

void bounding_rectangle(rectangle *a, rectangle *b, rectangle *result){
  result->x1=MIN(a->x1,b->x1);
  result->x2=MAX(a->x2,b->x2);
  result->y1=MIN(a->y1,b->y1);
  result->y2=MAX(a->y2,b->y2);
}

void inline shift_rect(int shift_x, int shift_y, rectangle *rect){
  rect->x1+=shift_x;rect->x2+=shift_x;
  rect->y1+=shift_y;rect->y2+=shift_y;
}

void adjust_rectangle(rectangle *a, rectangle *b, rectangle *result){
  /*Return bounding rectangle starting at (0,0). Modify rectangles a and b*/
  bounding_rectangle(a,b,result);
  int shift_x=-result->x1;
  int shift_y=-result->y1;
  shift_rect(shift_x,shift_y,result);
  shift_rect(shift_x,shift_y,a);
  shift_rect(shift_x,shift_y,b);
}  

void print_rectangle(rectangle *rect){
  printf("(%d,%d)--(%d,%d)",rect->x1,rect->y1,rect->x2,rect->y2);
}

double inline convertrgb(DATA32 pix){
  unsigned char r,g,b;
  b=pix & 0xff;
  g=(pix>>8) & 0xff;
  r=(pix>>16) & 0xff;
  return 1.0-(0.2126*r+0.7152*g+0.0722*b)/255.0;
}

DATA32 inline convertfromdouble(double pixvalue){
  unsigned char r,g,b;
  r=(int)(255.0*(1.0-pixvalue/50000));
  g=b=r;
  return (r<<16)|(g<<8)|b;
}

double redpart(DATA32 pix){
  unsigned char r;
  r=(pix>>16) & 0xff;
  return (double)r/255.0;
}

double greenpart(DATA32 pix){
  unsigned char r;
  r=(pix>>8) & 0xff;
  return (double)r/255.0;
}

double bluepart(DATA32 pix){
  unsigned char r;
  r=(pix) & 0xff;
  return (double)r/255.0;
}

fftw_plan plan, reverse_plan;
double *fftw_input;
fftw_complex *fftw_output;
double *fftw_reverse_result;

void initialize_fftw_engine(int size_x, int size_y){
  fftw_input=(double*)fftw_malloc(sizeof(double)*size_x*size_y);
  fftw_output=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*size_x*(size_y/2+1));
  fftw_reverse_result=(double*) fftw_malloc(sizeof(double)*size_x*size_y);
  if((!fftw_input)||(!fftw_output)||(!fftw_reverse_result)) die("Malloc failed");
  plan=fftw_plan_dft_r2c_2d(size_x,size_y,fftw_input,fftw_output,FFTW_MEASURE);
  reverse_plan=fftw_plan_dft_c2r_2d(size_x,size_y,fftw_output,fftw_reverse_result,FFTW_MEASURE);
  if(!plan) die("Initialisation of fft engine failed");
  if(!reverse_plan) die("Initialisation of fft engine failed");
}

void shutdown_fftw(){
  fftw_destroy_plan(plan);
  fftw_destroy_plan(reverse_plan);
  fftw_free(fftw_input);
  fftw_free(fftw_output);
  fftw_free(fftw_reverse_result);
}

void fftw_image_region(int size_x, int size_y, int offset_x, int offset_y, Imlib_Image image){
  /*Depends on the globals: plan, fftw_input, fftw_output (via plan) */
  /* The above mentioned globals should be preallocated and suitable to store the results */
  /* size_x and size_y parameters are actually redundant and should agree with plan and etc. */
  imlib_context_set_image(image);
  
  int width,height;
  width=imlib_image_get_width();
  height=imlib_image_get_height();
  if((size_x+offset_x>width)||(size_y+offset_y>height)) die("Transform region doesn't fit into image canvas");
  
  DATA32* image_bytes=imlib_image_get_data_for_reading_only();
  int x,y;
  for(y=0;y<size_y; y++){
   for(x=0; x<size_x; x++){
    fftw_input[y*size_x+x]=convertrgb(image_bytes[(y+offset_y)*width+x+offset_x]);
   }
  } 

  fftw_execute(plan);
}  
   
void save_image_region(int size_x, int size_y, int offset_x, int offset_y, Imlib_Image image,const char* filename){
  /*Debug function. */
  /* Stores region of image with given sizes and offsets into filename */
  /* Can be used to verify correctness of the offsets */
  imlib_context_set_image(image);
  
  int width,height;
  width=imlib_image_get_width();
  height=imlib_image_get_height();
  if((size_x+offset_x>width)||(size_y+offset_y>height)) die("Transform region doesn't fit into image canvas");
  
  DATA32* image_bytes=imlib_image_get_data_for_reading_only();
  DATA32* outbuffer=fftw_malloc(size_x*size_y*sizeof(DATA32));
  int x,y;
  for(y=0;y<size_y; y++){
   for(x=0; x<size_x; x++){
     outbuffer[y*size_x+x]=image_bytes[(y+offset_y)*width+x+offset_x];
   }
  }

  Imlib_Image out_image;
  out_image=imlib_create_image_using_data(size_x, size_y,outbuffer);
  imlib_context_set_image(out_image);
  imlib_image_set_format("png");
  imlib_save_image(filename); 
  imlib_free_image();
  fftw_free(outbuffer);
}

/*Functions for command-line arguments */
const char *progname;
void usage(){
  fprintf(stderr,"Usage %s [-s[taple] geom ] image1 image2 output_file\n ",progname);
  fprintf(stderr,"  where geom = [t|b|l|r]XxY, e.g. b200x200. The above examlple b200x200 means try to staple bottom of the first page to the top of the second\n");
  exit(-1);
}

int die(char *message)
{
  fprintf(stderr,"%s\n",message);
  exit(-1);
}

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



int main(int argc, const char **argv){
  int i;
  progname=argv[0];
  int debug_fft=0;
  Imlib_Image image1,image2;
  clock_t time_for_malloc, time_for_plan, time_for_load, time_for_array_access, time_for_fft;
  clock_t start_clock, end_clock;
  
  /*Parse arguments*/
  int size_x,size_y;
  char staple='b';
  const char *file[3]; //File names
  int filec=0; //Count of filenames
  char* geometry=NULL;
  for(i=1; i<argc; i++){
    const char *arg=argv[i];
    if(arg[0]=='-'){
      switch(arg[1]){
        case 's' :
          if(++i>=argc) usage(); //Increment i anyway
          geometry=(char*)argv[i];      //to consume next argument
          break;
        case 'v' :
          debug_fft=1;
          break;
        default:  
          usage();
      }
    }else{
      file[filec]=argv[i];
      filec++;
    }
  }  
  if(filec<3)
     usage();
  if(!geometry){
     size_x=SIZE_X; size_y=SIZE_Y;
     staple='b';
  }else{
    if(!parse_geometry(geometry, &size_x, &size_y, &staple)) usage();
  }   




  initialize_fftw_engine(size_x,size_y);
  double complex *storage_buffer=(double complex*)fftw_malloc(sizeof(double complex)*size_x*size_y);

  /*Load images*/
  image1=imlib_load_image(file[0]);
  if(!image1) die("Error loading first image");
  int width1,height1,width2,height2;
  imlib_context_set_image(image1);
  width1=imlib_image_get_width();
  height1=imlib_image_get_height();

  image2=imlib_load_image(file[1]);
  if(!image2) die("Error loading second image");
  imlib_context_set_image(image2);
  width2=imlib_image_get_width();
  height2=imlib_image_get_height();


  int off_x_1,off_y_1,off_x_2,off_y_2;
  switch(staple){
    case 'b' :
     off_x_1=(width1-size_x)/2;
     off_y_1=height1-size_y;
     off_x_2=(width2-size_x)/2;
     off_y_2=0;
     break;
    case 't' :
     off_x_1=(width1-size_x)/2;
     off_y_1=0;
     off_x_2=(width2-size_x)/2;
     off_y_2=height2-size_y;
     break;
    case 'r' :
     off_x_1=width1-size_x;
     off_y_1=(height1-size_y)/2;
     off_x_2=0;
     off_y_2=(height2-size_y)/2;
     break;
    case 'l' :
     off_x_1=0;
     off_y_1=(height1-size_y)/2;
     off_x_2=width2-size_x;
     off_y_2=(height2-size_y)/2; 
    default:
     off_x_1=0;
     off_x_2=0;
     off_y_1=0;
     off_y_2=0;
   } 

  fftw_image_region(size_x,size_y,off_x_1,off_y_1,image1);
  if(debug_fft)
    save_image_region(size_x,size_y,off_x_1,off_y_1,image1,"/tmp/test1.png");
  memcpy(storage_buffer,fftw_output,sizeof(double complex)*size_x*(size_y/2+1));

  //fftw_execute(reverse_plan);

  fftw_image_region(size_x,size_y,off_x_2,off_y_2,image2);
  if(debug_fft)
    save_image_region(size_x,size_y,off_x_2,off_y_2,image2,"/tmp/test2.png");
  for(i=0;i<size_x*(size_y/2+1);i++){
    /*See http://en.wikipedia.org/wiki/Phase_correlation */
    double complex tmp=conj(fftw_output[i])*storage_buffer[i];
    fftw_output[i]=tmp/cabs(tmp);
  }
  fftw_execute(reverse_plan);


  /*Find max of fft transformed correlation -- will indicate shift */ 
  double max=fftw_reverse_result[0];
  int imax=0;
  for(i=0; i<size_x*size_y;i++)
    if(fabs(fftw_reverse_result[i])>max){ imax=i; max=fabs(fftw_reverse_result[i]); }
  
  int shift_x, shift_y;
  shift_x=imax%size_x;
  shift_y=(imax-shift_x)/size_x;
  if(shift_x>size_x/2)
    shift_x=-(size_x-shift_x);
  if(shift_y>size_y/2)
    shift_y=-(size_y-shift_y);

  if(debug_fft)
    printf("Shift x: %d Shift y: %d\n",shift_x,shift_y);
  shift_y+=(off_y_1-off_y_2);
  shift_x+=(off_x_1-off_x_2);
  if(debug_fft)
    printf("Shift x: %d Shift y: %d\n",shift_x,shift_y);

  rectangle rect_im1, rect_im2, bbox;
  make_rect(0,0,width1,height1,&rect_im1);
  make_rect(0,0,width2,height2,&rect_im2);
  shift_rect(shift_x,shift_y,&rect_im2);
  if(debug_fft){
    printf("Before transform:\n");
    printf("Image 1: "); print_rectangle(&rect_im1); printf("\n");
    printf("Image 2: "); print_rectangle(&rect_im2); printf("\n");
  }
  adjust_rectangle(&rect_im1,&rect_im2,&bbox);
  if(debug_fft){
    printf("After transform:\n");
    printf("Image 1: "); print_rectangle(&rect_im1); printf("\n");
    printf("Image 2: "); print_rectangle(&rect_im2); printf("\n");
    printf("BBox: "); print_rectangle(&bbox); printf("\n");
  }

  
  Imlib_Image result_image;
  result_image=imlib_create_image(bbox.x2,bbox.y2);
  if(!result_image) die("Unable to create output image");
  imlib_context_set_image(result_image);
  imlib_image_set_format("png");
  imlib_context_set_color(255,255,255,255);
  imlib_image_fill_rectangle(0,0,bbox.x2,bbox.y2);
  imlib_blend_image_onto_image(image1,1,0,0,width1,height1,rect_im1.x1,rect_im1.y1,width1,height1);
  imlib_blend_image_onto_image(image2,1,0,0,width2,height2,rect_im2.x1,rect_im2.y1,width2,height2);
  imlib_save_image(file[2]);
  imlib_free_image();

  imlib_context_set_image(image1);
  imlib_free_image();
  imlib_context_set_image(image2);
  imlib_free_image();

  if(debug_fft){
    DATA32* outbuffer=fftw_malloc(sizeof(DATA32)*size_x*size_y);
    for(i=0; i<size_x*size_y; i++)
      outbuffer[i]=convertfromdouble(fabs(fftw_reverse_result[i]));

    Imlib_Image out_image;
    out_image=imlib_create_image_using_data(size_x, size_y,outbuffer);
    imlib_context_set_image(out_image);
    imlib_image_set_format("png");
    imlib_save_image("/tmp/tmp_fft.png"); 
    imlib_free_image();
    fftw_free(outbuffer);
  }


}  

