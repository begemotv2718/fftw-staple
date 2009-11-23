#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))


typedef struct{
  int x1,y1,x2,y2;
} rectangle;

void inline make_rect(int x1, int y1, int x2, int y2, rectangle *rect){
  rect->x1=x1;
  rect->y1=y1;
  rect->x2=x2;
  rect->y2=y2;
}

void inline copy_rect(rectangle *dest, rectangle *src){
  dest->x1=src->x1;
  dest->x2=src->x2;
  dest->y1=src->y1;
  dest->y2=src->y2;
}

void inline bounding_rectangle(rectangle *a, rectangle *b, rectangle *result){
  result->x1=MIN(a->x1,b->x1);
  result->x2=MAX(a->x2,b->x2);
  result->y1=MIN(a->y1,b->y1);
  result->y2=MAX(a->y2,b->y2);
}

void inline shift_rect(int shift_x, int shift_y, rectangle *rect){
  rect->x1+=shift_x;rect->x2+=shift_x;
  rect->y1+=shift_y;rect->y2+=shift_y;
}

int inline rect_width(rectangle rect){
  return rect.x2-rect.x1;
}

int inline rect_height(rectangle rect){
  return rect.y2-rect.y1;
}

void inline adjust_rectangle(rectangle *a, rectangle *b, rectangle *result){
  /*Return bounding rectangle starting at (0,0). Modify rectangles a and b*/
  bounding_rectangle(a,b,result);
  int shift_x=-result->x1;
  int shift_y=-result->y1;
  shift_rect(shift_x,shift_y,result);
  shift_rect(shift_x,shift_y,a);
  shift_rect(shift_x,shift_y,b);
}  

void inline print_rectangle(rectangle *rect){
  printf("(%d,%d)--(%d,%d)",rect->x1,rect->y1,rect->x2,rect->y2);
}
