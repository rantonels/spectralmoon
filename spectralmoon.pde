int sz = 640;

float maxl = sz * sqrt(2);

int iters = 10;

float ts = 0.3;

float [][] walls = {
  {0,0,sz,0},
  {sz,0,sz,sz},
  {sz,sz,0,sz},
  {0,sz,0,0},

  {0,0,0,0},
  {0,0,0,0},
  {0,0,0,0},
  {0,0,0,0},
  {0,0,0,0},
  {0,0,0,0},
  {0,0,0,0},
  {0,0,0,0},
  {0,0,0,0}
};

int WALL = 0;
int MIRROR = 1;
int GLASS = 2;

int type [] = {
  0,
  0,
  0,
  0,
  
  2,
  2,
  2,
  2,
  2,
  2,
  2,
  2,
  2
};

int nwalls = walls.length;

int segments = 70;
int nrays = 70;
PGraphics [] pg = new PGraphics[segments];  

void triangle(float theta)
{
 for (int i=0; i<9;i++)
  {
    float angle = theta+ -PI/18  + 0.01 + TWO_PI / 9. * i;
    float angle2 = theta + -PI/18 + 0.01 + TWO_PI / 9. * ((i+1)%9);
    float r = sz * ts;
    walls[4+i][0] = sz/2 + r*cos(angle);
    walls[4+i][1] = sz/2 - r*sin(angle);
    walls[4+i][2] = sz/2 + r*cos(angle2);
    walls[4+i][3] = sz/2 - r*sin(angle2);
  }
}

void setup()
{
  
  
  size(sz,sz,P2D);
  for (int i =0; i<segments; i++)
    pg[i] = createGraphics(sz,sz,P2D);
  colorMode(RGB,255);
 frameRate(24);
  
 
};

float n(float lambda)
{
  //return 1.7280 + 13420/(lambda*lambda);
  return 1.435 + 30000/(lambda*lambda);

}

//from http://stackoverflow.com/questions/3407942/rgb-values-of-visible-spectrum

color spectral_color(float l,float alpha) // RGB <0,1> <- lambda l <400,700> [nm]
    {
    float t;  
    float r=0.0; 
    float g=0.0; 
    float b=0.0;
    
         if ((l>=400.0)&(l<410.0)) { t=(l-400.0)/(410.0-400.0); r=    +(0.33*t)-(0.20*t*t); }
    else if ((l>=410.0)&(l<475.0)) { t=(l-410.0)/(475.0-410.0); r=0.14         -(0.13*t*t); }
    else if ((l>=545.0)&(l<595.0)) { t=(l-545.0)/(595.0-545.0); r=    +(1.98*t)-(     t*t); }
    else if ((l>=595.0)&(l<650.0)) { t=(l-595.0)/(650.0-595.0); r=0.98+(0.06*t)-(0.40*t*t); }
    else if ((l>=650.0)&(l<700.0)) { t=(l-650.0)/(700.0-650.0); r=0.65-(0.84*t)+(0.20*t*t); }
         if ((l>=415.0)&(l<475.0)) { t=(l-415.0)/(475.0-415.0); g=             +(0.80*t*t); }
    else if ((l>=475.0)&(l<590.0)) { t=(l-475.0)/(590.0-475.0); g=0.8 +(0.76*t)-(0.80*t*t); }
    else if ((l>=585.0)&(l<639.0)) { t=(l-585.0)/(639.0-585.0); g=0.84-(0.84*t)           ; }
         if ((l>=400.0)&(l<475.0)) { t=(l-400.0)/(475.0-400.0); b=    +(2.20*t)-(1.50*t*t); }
    else if ((l>=475.0)&(l<560.0)) { t=(l-475.0)/(560.0-475.0); b=0.7 -(     t)+(0.30*t*t); };
    
    int R = int(255.*r*alpha);
    int G = int(255.*g*alpha);
    int B = int(255.*b*alpha);
    
    return color(R,G,B);
  }
 
PVector intersection(PVector p1, PVector p2, PVector p3, PVector p4) {
// Store the values for fast access and easy
// equations-to-code conversion
float x1 = p1.x, x2 = p2.x, x3 = p3.x, x4 = p4.x;
float y1 = p1.y, y2 = p2.y, y3 = p3.y, y4 = p4.y;
 
float d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
// If d is zero, there is no intersection
if (d == 0) return null;
 
// Get the x and y
float pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
float x = ( pre * (x3 - x4) - (x1 - x2) * post ) / d;
float y = ( pre * (y3 - y4) - (y1 - y2) * post ) / d;
 
// Check if the x and y coordinates are within both lines
if ( x < min(x1, x2) | x > max(x1, x2) |
x < min(x3, x4) | x > max(x3, x4) ) return null;
if ( y < min(y1, y2) | y > max(y1, y2) |
y < min(y3, y4) | y > max(y3, y4) ) return null;
 
// Return the point of intersection
PVector ret = new PVector(x,y);

return ret;
} 
 
float clip(float t)
{
  return min(1,max(0,t));
};

void trace_ray(float alpha, PGraphics pgp)
{
  pgp.blendMode(ADD);
  PVector pos = new PVector(0.5,330*sz/float(640)+0.0005*randomGaussian()*sz);
  PVector ray = new PVector(1,-0.25*(2.4+sin(frameCount/60.*TWO_PI))-0.002*randomGaussian());
  

  
  ray.normalize();
  float wl = 400 + (700-400)*random(1000)/1000.;
  //float wl = 600;
  
  for (int it=0; it < iters; it++)
  {
    //wl+=50;
    PVector inters = PVector.add(pos,PVector.mult(ray,maxl )); 
    PVector finalp = inters.get();
    
    float distance = maxl;

    int intwall = -1;    
    
    for (int j=0; j<nwalls; j++)
    {
      PVector w1 = new PVector(walls[j][0],walls[j][1]);
      PVector w2 = new PVector(walls[j][2],walls[j][3]);
      PVector jint = intersection(pos,finalp,w1,w2);
      
      if (jint != null)
      {
        float tmpdistance = PVector.dist(pos,jint);
        if(tmpdistance < distance)
        {
          //good intersection
          
          distance = tmpdistance;
          inters = jint;
          intwall = j;
        };
      };
      
    };
    
   
    
    pgp.strokeWeight(1);
    pgp.stroke(spectral_color(wl,alpha));
    pgp.line(pos.x,pos.y,inters.x,inters.y);
    
    
    
    if (intwall == -1)
      break;
      
    //fill(255,0,0,alpha*0.3);
    //ellipse(inters.x,inters.y,10*it,10*it);
      
    if(type[intwall] == 0) //no intersection / solid wall
      break;
    
    PVector w1 = new PVector(walls[intwall][0],walls[intwall][1]);
    PVector w2 = new PVector(walls[intwall][2],walls[intwall][3]);
    PVector lon = PVector.sub(w2,w1);
    lon.normalize();
    
    //(1 0) -> (0 -1)
    //(0 1) -> (1  0)

    PVector normal = new PVector(-lon.y,lon.x);
    
    if (type[intwall] == MIRROR) //mirror reflection
    {
      
      ray = PVector.sub(ray , PVector.mult(normal,2*PVector.dot(normal,ray)));
      
      ray.normalize();
      
      pos = PVector.add(inters.get(),PVector.mult(ray,1));
      
    }
    
    if (type[intwall] == GLASS)
    {
      ray.normalize();
      PVector nnorm = normal.get();
      nnorm.normalize();
      boolean isswap = false;
      if(PVector.dot(normal,ray) > 0)
      {
        nnorm = PVector.mult(nnorm,-1);
        isswap = true;
        
      }
      
      float n = n(wl);
      float nn = n;
      float sgn = 1.;
      if (isswap)
      {
          nn = 1./nn;
          sgn = -1.;
      }
     float costhetai = abs(PVector.dot(ray,nnorm));
     
     float sinthetatsq = pow(nn,-2) * (1-costhetai*costhetai); 
     
     float reflectance;
     if(sinthetatsq>1)
     {
      //print("total reflection");
      reflectance = 1.;
     }
     else
     {
     float costhetat = sqrt(1-sinthetatsq);
    

     float Rs_sqrt = (costhetai - nn*costhetat)/(costhetai + nn*costhetat);
     
     float Rp_sqrt = (costhetat - nn*costhetai)/(costhetat + nn*costhetai);
     
     float Rs = Rs_sqrt * Rs_sqrt;
     float Rp = Rp_sqrt * Rp_sqrt;

     reflectance = .5*Rp*Rs; 
     //print(reflectance);   
     }
      
      if(random(10000)/10000. < reflectance) //reflection
      {
        ray = PVector.sub(ray , PVector.mult(nnorm,2*PVector.dot(nnorm,ray))); 
        ray.normalize();
        pos = PVector.add(inters.get(),PVector.mult(ray,1));
        
      }
      else      //refraction
      {
        
        
        ray = PVector.sub(
              PVector.mult(ray,1/nn),
              PVector.mult(nnorm,(-1/nn*costhetai + sqrt(1-sinthetatsq)))
              );
              
        ray.normalize();
        
        pos = PVector.add(inters.get(),PVector.mult(ray,2));
      };
      
    };
    
    
    
    
    
  };
  print("\n");
}
    
void draw_debug()
{
  stroke(255);
  strokeWeight(4);
  for(int j=0; j<nwalls; j++)
  {
    line(walls[j][0],walls[j][1],walls[j][2],walls[j][3]);
//    PVector w1 = new PVector(walls[j][0],walls[j][1]);
//    PVector w2 = new PVector(walls[j][2],walls[j][3]);
//    PVector lon = PVector.sub(w2,w1);
//    lon.normalize();
//    
//    //(1 0) -> (0 -1)
//    //(0 1) -> (1  0)
//
//    PVector normal = new PVector(-lon.y,+lon.x);
//    normal.normalize();
//    
//    line(walls[j][0],walls[j][1],walls[j][0]+20*normal.x,walls[j][1]+20*normal.y);
//    
  }  
}
   
float smootherstep(float x)
{

   return x*x*x*(x*(x*6-15)+10);
}
   
float cut = 0.5;
   
void draw()
{
 float progress = 1/60.*(frameCount%60);
 float theta;
 if (progress < cut)
   theta = 0;
 else
   theta = TWO_PI/3. * smootherstep((progress-cut)/(1-cut));
   //theta = 0;
   
 triangle(theta);
  
 background(0);
 for (int i = 0; i<segments; i++)
 {
   pg[i].beginDraw();
   pg[i].blendMode(BLEND);
   pg[i].background(0);

   pg[i].blendMode(ADD);
   for (int j = 0; j<nrays; j++)
   {
     trace_ray(1/float(nrays),pg[i]);
   }
   
   pg[i].blendMode(BLEND);
   pg[i].fill(0,0,0,float(segments-1)/float(segments));
   pg[i].rect(0,0,sz,sz);
   
   pg[i].endDraw();
   
   blendMode(ADD);
   image(pg[i],0,0);
 };

 draw_debug();
 
 //if (frameCount < 60)
 //  saveFrame("nonagon-##.png");
};
