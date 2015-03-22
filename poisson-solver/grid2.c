#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#endif
#define GL_GLEXT_PROTOTYPES
#include <GL/glut.h>

int main()
{
  //Create a 1D image - for this example it's just a red line
  int stripeImageWidth = 32;
  GLubyte   stripeImage[3*stripeImageWidth];
  for (int j = 0; j < stripeImageWidth; j++) {
      stripeImage[3*j] = j*255/32; // use a gradient instead of a line
      stripeImage[3*j+1] = 255;
      stripeImage[3*j+2] = 255;
  }
  
  GLuint texID;
  glGenTextures(1, &texID);
  glBindTexture(GL_TEXTURE_1D, texID);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage1D(GL_TEXTURE_1D, 0, 3, stripeImageWidth, 0, GL_RGB, GL_UNSIGNED_BYTE, stripeImage);
  
  // We want the texture to wrap, so that values outside the range [0, 1] 
  // are mapped into a gradient sawtooth
  glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

  // The texture coordinate comes from the data, it it not
  // generated from the vertex position!!!
  glDisable( GL_TEXTURE_GEN_S ); 
  glDisable(GL_TEXTURE_2D);
  glEnable( GL_TEXTURE_1D );
  glBindTexture(GL_TEXTURE_1D, texID);
  
  float hist2D[bins_x][bins_y] = {NaN, NaN, ...}
  
  glBegin(GL_QUADS);
  for(int y=0; y<grid_height; y+=2) for(int x=0; x<grid_width; x+=2) {
      glTexCoord1f(hist2D[x  ][y  ]]); glVertex2i(x  ,y);
      glTexCoord1f(hist2D[x+1][y  ]]); glVertex2i(x+1,y);
      glTexCoord1f(hist2D[x+1][y+1]]); glVertex2i(x+1,y+1);
      glTexCoord1f(hist2D[x  ][y+1]]); glVertex2i(x  ,y+1);
  }
  glEnd();
  
  return 0;
}