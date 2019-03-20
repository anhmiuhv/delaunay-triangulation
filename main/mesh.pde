import java.util.*;
// TRIANGLE MESH
class MESH {
    // VERTICES
    int nv=0, maxnv = 1000;  
    pt[] G = new pt [maxnv];
    //QUADS
    int nq = 0, maxnq = maxnv;    
    // TRIANGLES 
    int nt = 0, maxnt = maxnv*2;                           
    boolean[] isInterior = new boolean[maxnv];                                      
    // CORNERS 
    int c=0;    // current corner
    int C=0; //current quad corner
    int nc = 0;
    int nC = 0;
    int[] V = new int [3*maxnt];   
    int[] O = new int [3*maxnt];
    Integer[] S = new Integer [4*maxnq];
    Integer[] M = new Integer [3*maxnt];
    Integer[] P = new Integer [3*maxnt];
    // current corner that can be edited with keys
  MESH() {for (int i=0; i<maxnv; i++) G[i]=new pt();};
  void reset() {nv=0; nq=0; nt=0; nc=0; nC=0;}                                                  // removes all vertices and triangles
  void loadVertices(pt[] P, int n) {nv=0; for (int i=0; i<n; i++) addVertex(P[i]);}
  void writeVerticesTo(pts P) {for (int i=0; i<nv; i++) P.G[i].setTo(G[i]);}
  void addVertex(pt P) { G[nv++].setTo(P); }                                             // adds a vertex to vertex table G
  void addTriangle(int i, int j, int k) {V[nc++]=i; V[nc++]=j; V[nc++]=k; nt=nc/3; }     // adds triangle (i,j,k) to V table
  void addQuad(int i, int j, int k, int l) {V[nC++] = i; V[nC++]=j; V[nC++]=k; V[nC++]=l; nq=nC/4;} // adds Quad   (i,j,k,l) to V table
  
  //Quad Corner Operations
  int C (int V) {return 4*V;}
  int V (int C) {
    if (C%4 == 0 && C < 4*nv)
      return C/4;
    return V(S(C));
  }
  int Q(int C){return C/4;}
  int N(int C){return 4*Q(C) + ((C+1)%4);}
  int S(int C){return S[C];}
  
  //Corner Conversions
  int C(int Q, int i) {return 4*Q + i;}
  int c(int Q, int i) {return 8*Q + i;}
  int Qc(int c){return c/8;}
  int ConvertC(int c){return C(Qc(c),(c+ ((c/2)&2))%4);}
  boolean isQuad(int Q){return S[4*Q+3]!=null;}
  int Convertc(int C){
    if (isQuad(C))
      return c(Q(C),(C%4)+(C&2));
    return c(Q(C),C%4);
  }
  //Triangle Meshes
  int c(int v) {return 8*v;}
  int ct(int t) {return 4*t;}
  int v(int c){return V(ConvertC(c));}
  int t(int c){return c/4;}
  int n(int c){
    if (c%4 == 2)
      return c-2;
    return c+1;
  }
  int s(int c){
    if (c%8 == 0 && isQuad(Qc(c)))
      return c(Qc(c),6);
    if (c%8 == 4)
      return c(Qc(c),2);
    return Convertc(S(ConvertC(c)));
  }
  
  
  // CORNER OPERATORS
  //int t (int c) {int r=int(c/3); return(r);}                   // triangle of corner c
  //int n (int c) {int r=3*int(c/3)+(c+1)%3; return(r);}         // next corner
  int p (int c) {return n(n(c));}//{int r=3*int(c/3)+(c+2)%3; return(r);}         // previous corner
  pt g (int c) {return G[V[c]];}                             // shortcut to get the point where the vertex v(c) of corner c is located

  boolean nb(int c) {return(o(c)!=c);};  // not a border corner
  boolean bord(int c) {return(o(c)==c);};  // not a border corner

  pt cg(int c) {return P(0.6,g(c),0.2,g(p(c)),0.2,g(n(c)));}   // computes offset location of point at corner c

  // CORNER ACTIONS CURRENT CORNER c
  void next() {c=n(c);}
  void previous() {c=p(c);}
  void opposite() {c=o(c);}
  void left() {c=l(c);}
  void right() {c=r(c);}
  void swing() {c=s(c);} 
  void unswing() {c=u(c);} 
  void printCorner() {println("c = "+c);}
  
  

  // DISPLAY
  void showCurrentCorner(float r) { if(bord(c)) fill(red); else fill(dgreen); show(cg(c),r); };   // renders corner c as small ball
  void showEdge(int c) {beam( g(p(c)),g(n(c)),rt ); };  // draws edge of t(c) opposite to corner c
  void showVertices(float r) // shows all vertices green inside, red outside
    {
    for (int v=0; v<nv; v++) 
      {
      if(isInterior[v]) fill(green); else fill(red);
      show(G[v],r);
      }
    }                          
  void showInteriorVertices(float r) {for (int v=0; v<nv; v++) if(isInterior[v]) show(G[v],r); }                          // shows all vertices as dots
  void showTriangles() { for (int c=0; c<nc; c+=3) show(g(c), g(c+1), g(c+2)); }         // draws all triangles (edges, or filled)
  void showEdges() {for (int i=0; i<nc; i++) showEdge(i); };         // draws all edges of mesh twice

  void triangulate()      // performs Delaunay triangulation using a quartic algorithm
   {
   c=0;                   // to reset current corner
   // **01 implement it
   for (int i=0; i<nv; i++){
     for (int j=i; j<nv; j++){
       for (int k=j; k<nv; k++){
         if(i != j && i != k && j!= k){
           if(!isFlatterThan(G[i],G[j],G[k],10.0)){
             pt C = CircumCenter(G[i],G[j],G[k]);
             boolean good = true;
             for (int m=0; m<nv; m++){
               if(m!=i && m!=j && m!=k){
                 if(d(C,G[m]) < d(C,G[i])){
                   good = false;                      
                 }
               }
             }
             if(good){
               if(ccw(G[i],G[j],G[k])){
                 addTriangle(i,j,k);
               }
               else{
                 addTriangle(i,k,j);
               }
             }
           }
         }
       }
     }
   }
   }  

  void MatchAndPair(int c){
    M[v(c)] = t(c);
    M[v(n(c))] = t(s(n(c)));
    M[v(p(c))] = t(s(p(c)));
    P[t(s(n(c)))] = t(s(n(c)));
    P[t(s(p(c)))] = t(s(p(c)));
    boolean [] T = new boolean[nt];
    T[t(c)] = true;
    c = l(c);
    Stack<Integer> stack = new Stack<Integer>();
    stack.push(null);
    while (!stack.empty()){
      T[t(c)] = true;
      if (M[v(c)] == null){
        M[v(c)] = t(c);
        if (P[t(r(c))] == null && M[v(r(c))] != null){
          P[t(c)] = t(r(c));
          P[t(r(c))] = t(c);
        }
        else if (P[t(l(c))] == null){
          P[t(c)] = t(l(c));
          P[t(l(c))] = t(c);
        }
        c = r(c);
      }
      else{
        if (T[t(l(c))] == true){
          if(T[t(r(c))] == true)
            c = stack.pop();
          else
            c = r(c);
        }
        else {
          if(T[t(r(c))] == true)
            c = l(c);
          else {
            stack.push(l(c));
            c = r(c);
          }
        }
      }
    }
  }
  
  
   
  void computeO() // **02 implement it 
    {                                          
    // **02 implement it
    for (int i=0; i<nc; i++){
      boolean hasOpposite = false;
      for (int j=0; j<nc; j++){
        if(i != j){
          if(g(n(i)) == g(p(j)) && g(p(i)) == g(n(j))){
            O[i] = j;
            hasOpposite = true;
          }
        }
      }
      if(!hasOpposite){
        O[i] = i;
      }
    }
    } 
    
  void showBorderEdges()  // draws all border edges of mesh
    {
    // **02 implement;
    for (int i=0; i<nc; i++){
      if(bord(i)){
        showEdge(i);
      }
    }
    }
  
  int countBorders()
    {
      int count = 0;
      for (int i=0; i<nc; i++){
        if(bord(i)){
          count++;
        }
      }
      return count;
    }

  void showNonBorderEdges() // draws all non-border edges of mesh
    {
    // **02 implement
    for (int i=0; i<nc; i++){
      if(nb(i)){
        showEdge(i);
      }
    }
    }        
    
  void classifyVertices() 
    { 
    // **03 implement it
    for (int i=0; i<nv; i++){
      boolean Interior = true;
      for (int j=0; j<nc; j++){
        if(V[j] == i){
          if(bord(n(j)) || bord(p(j))){
            Interior = false;
          }
        } 
      }
      isInterior[i] = Interior;
    }
    }  
    
  void smoothenInterior() { // even interior vertiex locations
    pt[] Gn = new pt[nv];
    // **04 implement it 
    for (int v=0; v<nv; v++){
      int count = 0;
      for (int c=0; c<nc; c++){
        if(V[c] == v){
          if(Gn[v] == null){
            Gn[v] = g(n(c));
            Gn[v] = A(Gn[v],g(p(c)));
            count += 2;
          }
          else{
            Gn[v] = A(Gn[v],g(n(c)));
            Gn[v] = A(Gn[v],g(p(c)));
            count += 2;            
          }
        }
      }
      Gn[v].x = Gn[v].x / count;
      Gn[v].y = Gn[v].y / count;
      Gn[v].z = Gn[v].z / count;
    }
    for (int v=0; v<nv; v++){ 
      if(isInterior[v]){        
        G[v].translateTowards(.1,Gn[v]);
      }
    }
    }


   // **05 implement corner operators in Mesh
  //int v (int c) {return V[c];}                                // vertex of c
  int o (int c) {return p(s(p(c)));}                                // opposite corner
  int l (int c) {return o(n(c));}                             // left
  //int s (int c) {return n(l(c));}                             // left
  int u (int c) {return p(r(c));}                             // left
  int r (int c) {return o(p(c));}                             // right
  
  void showOpposites(){
    for (int c=0; c<nc; c++){
      if(nb(c)){
        drawParabolaInHat(cg(c), P(g(n(c)),g(p(c))), cg(o(c)), 2);
      }
    }
  }

  void showVoronoiEdges() // draws Voronoi edges on the boundary of Voroni cells of interior vertices
    { 
    // **06 implement it
    for (int v=0; v<nv; v++){
      if(isInterior[v]){
        for (int c=0; c<nc; c++){
          if(v(c) == v){
            int curr = c;
            int next;            
            do{
              next = s(curr);
              show(triCircumcenter(curr),triCircumcenter(next));
              curr = next;
            }while(curr != c);
            break;
          }
        }
      }
    }
    }               

  void showArcs() // draws arcs of quadratic B-spline of Voronoi boundary loops of interior vertices
    { 
    // **06 implement it
    for (int v=0; v<nv; v++){
      if(isInterior[v]){
        for (int c=0; c<nc; c++){
          if(v(c) == v){
            int curr = c;
            pt C1,C2;
            do{
              C1 = P(triCircumcenter(curr),triCircumcenter(s(curr)));
              C2 = P(triCircumcenter(s(curr)),triCircumcenter(s(s(curr))));
              drawParabolaInHat(C1, triCircumcenter(s(curr)), C2, 4);
              curr = s(curr);
            }while(curr != c);
            break;
          }
        }
      }
    }
    }               // draws arcs in triangles
  
  void drawVoronoiFaceOfInteriorVertices(){
    float dc = 1./(nv-1);
    for (int v=0; v<nv; v++){fill(dc*255*v,dc*255*(nv-v),200);drawVoronoiFaceOfInteriorVertices(v);}
  }
  void drawVoronoiFaceOfInteriorVertices(int v){
    if(isInterior[v]){
        for (int c=0; c<nc; c++){
          if(v(c) == v){
            int curr = c;
            pt C;
            beginShape();
            do{              
              C = triCircumcenter(curr);
              vertex(C.x,C.y,C.z);
              curr = s(curr);
            }while(curr != c);
            endShape(CLOSE);
            break;
          }
        }
      }
  }
 
  pt triCenter(int c) {return P(g(c),g(n(c)),g(p(c))); }  // returns center of mass of triangle of corner c
  pt triCircumcenter(int c) {return CircumCenter(g(c),g(n(c)),g(p(c))); }  // returns circumcenter of triangle of corner c


  } // end of MESH
