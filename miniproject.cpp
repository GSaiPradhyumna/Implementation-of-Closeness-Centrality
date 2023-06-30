#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <stack>

//#include <bits/stdc++.h>

using namespace std;

int inf=10e9;

class Edge {
    public :
    int from;
    int to;
};

struct node
{
    bool isartiPoint;
    int artiPoint;
    vector<int> vertices;
    node* parent;
    vector<node*> children;
    double DH;
};

class Graph 
{
    int v;// no of vertices
    vector <double> CC;       // Closeness centrality values of nodes.
    vector <int> deg;         //degreees of each vertex
    vector <int> arti;        // 1 if articulation point, esle 0.
    vector <int> inL;         // 1 if affected, else zero
    vector <int> inR;         // 1 if in R, else zero
    //vector <int> delH;        // values of delta of Hv 
    vector<vector<int> > d;   // Shortest distance matrix
    vector<vector<int> > prevd; 
    vector <double> D;           // D'[v]
    vector <double> prevD;       // D[v]
    vector<vector<int> > Adj;  // Adjacency matrix of the graph.
    node * rootBCT;

    vector <int> depth; // Depth of each vertex during DFS
    vector <int> lowLink; // Low-link value of each vertex during DFS
    vector <bool> visited; // Array to track visited vertices
    stack<pair <int, int> > st; // Stack to store the visited edges
    unordered_map<int,bool> articularPoints;// bool 1 for arti vertex else 0
    
    unordered_map<node* ,vector<int> > artiInbcc;
    unordered_map<int,vector<node*> > BCCofatri;
    unordered_map<int,node* > BCCofVert;
    //unordered_map<node*,int > visBCC; 
    unordered_map<node*,int > Delta; // map from BCC to its Delta value
    unordered_map<node*,int > bccno; // given a BCC node , gives its serial number
    unordered_map<int,node* > ap; // returns node* of an artiPoint
    vector<unordered_set<int> > bcc; // BCC content vertices
    vector<node*> BCC; // Only BCCs
    vector<vector<int > >delta ;// actual size of array is bcc.size(),v; 

    //map <pair<int, int>, int> dd;


    public:
    Graph (int vertNum) 
    {
        v=vertNum;

        CC.resize(v, 0);
        D.resize(v, 0);
        prevD.resize(v, 0);
        deg.resize(v, 0);
        inL.resize(v, 0);
        //delH.resize(v, 0);
        inR.resize(v, 0);
        arti.resize(v, 0);
        BCC.resize(v,NULL);

        depth.resize(v, 0);
        lowLink.resize(v, 0);
        visited.resize(v, 0);

        d.resize(v, vector<int>(v, inf));
        prevd.resize(v, vector<int>(v, inf));
        delta.resize(v,vector<int>(v,0));
        Adj.resize(v, vector<int>(v, inf));
        for (int i = 0; i < v; i++)
        {
            CC[i]=0;
            d[i][i] = 0;
            Adj[i][i] = 0;
            for (int j = 0; j < v; j++)
                Adj[i][j] = 0;
        }
    }    

    bool isR3(int vert)
    {
        if(deg[vert]!=3)
            return 0;
        int ee[3];
        int tt=0;
        for(int i=0;i<v;i++) 
        {
            if(Adj[vert][i]==1) 
                ee[tt++]=i;
        }
        for(int i=0;i<2;i++)
        {
            for(int j=i+1;j<3;j++)
            {
                if(Adj[ee[i]][ee[j]]==0)
                    return 0;
            }
        }
        return 1;
    }

    bool isR4(int vert)
    {
        if(deg[vert]!=4)
            return 0;
        int ee[4];
        int tt=0;
        for(int i=0;i<v;i++) 
        {
            if(Adj[vert][i]==1) 
                ee[tt++]=i;
        }
        int c=0;
        for(int i=0;i<4;i++)
        {
            c=0;
            for(int j=0;j<4;j++)
            {
                if(Adj[ee[i]][ee[j]]==1)
                    c++;
            }
            if(c<2)
                return 0;
        }
        return 1;
    }

    bool isChainNode(int vert)
    {
        if(deg[vert]!=2)
            return 0;
        int viss[v];
        int nb; //neighbour
        for(int i=0;i<v;i++)
        {
            viss[i]=0;
            if (Adj[vert][i]==1)
                nb=i;
        }
        viss[vert]=1;// visited
        int cc=0;
        while(true)
        {
            cc++;
            for(int i=0;i<v;i++)
            {
                if(deg[nb]!=2)
                    return 1;
                if(cc>1 && Adj[nb][i]==1 && i==vert) //chain occured
                    return 0;
                if(Adj[nb][i]==1 && viss[i]==0)
                {
                    viss[i]==1;
                    nb=i;
                    break;
                }
            }
        }
        return 1;
    }

    void findDegrees() {
        int c=0;
        for(int i=0;i<v;i++) 
        {
            c=0;
            for(int j=0;j<v;j++) 
            {
                if(Adj[i][j]==1)
                c++;
            }
            deg[i]=c;
        }
    }

    int BFS_Req(int vert, vector < int > SS)
    {
        int inSS[v];
        bool vis[v];
        for(int i=0;i<v;i++)
        {
            inSS[i]=0;
            vis[i]=false;
        }
        
        cout<<"SS : ";
        for(int i=0;i<SS.size();i++ )
        {
            cout<<SS[i]<<' ';
            inSS[SS[i]]=1;
        }
        cout<<endl;
        
        vis[vert]=true;
        
        
        queue <int>qq1;
        queue <int>qq2;
        qq1.push(vert);
        int u,dep=0,cc=0;
        int res=0;
        while(!qq1.empty())
        {
            dep++;
            u=qq1.front();
            qq1.pop();
            for(int i=0;i<v;i++)
            {
                if(Adj[u][i]==1 && vis[i]==false && inSS[i]==1)
                {
                    qq2.push(i);
                    vis[i]=true;
                }
            }
            res+=(dep*qq2.size());
            qq1=qq2;
            cout<<"q2 : ";
            while(!qq2.empty()){
                cout<<qq2.front()<<" ";
                qq2.pop();
            }
            cout<<endl;
                
        }
        return res;
    }
    
    void tarjanBCC(int vert, int parent) 
    {
        static int currentDepth = 0;
        depth[vert] = lowLink[vert] = ++currentDepth;
        int childCount = 0;
        visited[vert] = true;

        for (int i = 0; i < v; ++i) 
        {
            if(Adj[vert][i]==1)
            {
                int u=i;
                if (!visited[u]) {
                    //if not vissited to u then 
                    ++childCount;
                    st.push(make_pair(vert,u ));
                    tarjanBCC(u, vert);

                    lowLink[vert] = min(lowLink[vert], lowLink[u]);

                    if ((parent == -1 && childCount > 1) ) 
                    {
                        unordered_set<int> x;
                        while (st.top() != make_pair(vert, u)) {
                            x.insert(st.top().first);
                            x.insert(st.top().second);
                            cout << st.top().first << "-" << st.top().second << " ";
                            st.pop();
                        }
                        x.insert(st.top().first);
                        x.insert(st.top().second);
                        cout << st.top().first << "-" << st.top().second << endl;
                        st.pop();
                        articularPoints[vert]=1;
                        bcc.push_back(x);
                    }
                    else if ((parent == -1 && childCount == 1) ) 
                    {
                        unordered_set<int> x;
                        while (st.top() != make_pair(vert, u)) {
                            x.insert(st.top().first);
                            x.insert(st.top().second);
                            cout << st.top().first << "-" << st.top().second << " ";
                            st.pop();
                        }
                        x.insert(st.top().first);
                        x.insert(st.top().second);
                        cout << st.top().first << "-" << st.top().second << endl;
                        st.pop();
                        bcc.push_back(x);
                    }
                    else if((parent != -1 && lowLink[u] >= depth[vert]))
                    {
                        unordered_set<int> x;
                        while (st.top() != make_pair(vert, u)) {
                            x.insert(st.top().first);
                            x.insert(st.top().second);
                            cout << st.top().first << "-" << st.top().second << " ";
                            st.pop();
                        }
                        x.insert(st.top().first);
                        x.insert(st.top().second);
                        cout << st.top().first << "-" << st.top().second << endl;
                        st.pop();
                        articularPoints[vert]=1;
                        bcc.push_back(x);
                    }
                } else if (u != parent && depth[u] < lowLink[vert]) {
                    lowLink[vert] = min(lowLink[vert], depth[u]);
                    st.push(make_pair(vert,u ));
                }
            }
        }
    }

    void findBCC() 
    {
        // Input the adjacency matrix
        // Assuming the graph has n vertices (0 to n-1)
        // adjMatrix[i][j] = 1 if there is an edge between vertices i and j, 0 otherwise
        //int adjMatrix[MAX][MAX];



        // Read the number of vertices
        for (int i = 0; i < v; ++i) {
            depth[i] = lowLink[i] = visited[i] = 0;
            //adj[i].clear();
        }

        // Apply Tarjan's algorithm to find the BCCs
        cout << "Biconnected Components (BCCs):" << endl;
        for (int i = 0; i < v; ++i) {
            if (!visited[i]) {
                //cout<<" i = "<<i<<endl;
                tarjanBCC(i, -1);
            }
        }
        cout<<"HelloBRO"<<endl;
        return;
        cout<<bcc.size()<<"size"<<endl;
        cout<<"artipoints are ";
        for(int i=0;i<v;i++)
        {
            if(articularPoints[i]==1)
            {
                cout<<i<<" ";
            }
        }
        cout<<endl;
    }

    node* BCT()
    {
        // vector<node*> BCC;
        cout<<"pqrst"<<endl;
        // nodes corres to bcc are formed itr
        int b=0;
        cout<<"bcc.size() = "<<bcc.size()<<endl;

        for(int i=0;i<bcc.size();i++)
        {
            //cout<<"i = "<<i<<endl;
            b++;
            vector<int> v;
            struct node* temp = new node;
            temp->isartiPoint=false;
            temp->artiPoint=-1;
            temp->parent=NULL;
            unordered_set<int>::iterator itr;
            for( itr=bcc[i].begin();itr!=bcc[i].end();itr++)
            {
                int j=*itr;
                if(articularPoints[j]==1)
                {
                    BCCofatri[j].push_back(temp);
                    v.push_back(j);
                }
                else 
                {
                    BCCofVert[j]=temp;
                }
                (temp->vertices).push_back(j);
            }
            artiInbcc[temp]=v;
            BCC.push_back(temp);
        }
        //nodes corres to arti aare formed 
        int a=0;
        for(int i=0;i<v;i++)
        {
            if(articularPoints[i]==1)
            {
                struct node* temp = new node;
                temp->isartiPoint=true;
                temp->artiPoint=i;
                temp->parent=NULL;
                ap[i]=temp;
                a++;
            }
        }

        // combining bcc nodes and arti nodes
        queue<node*> que;
        que.push(BCC[0]);
        node* root=BCC[0];
        while(!que.empty())
        {
            node* t=que.front();
            que.pop();
            if(t->isartiPoint)
            {
                //articu point
                vector<node*> vec=BCCofatri[t->artiPoint];
                for(int i=0;i<vec.size();i++)
                {
                    if((vec[i]->children).size()==0)
                    {
                        (t->children).push_back(vec[i]);
                        vec[i]->parent=t;
                        que.push(vec[i]);
                    }
                }
            }
            else 
            {
                vector<int> vec=artiInbcc[t];
                for(int i=0;i<vec.size();i++)
                {
                    if(ap[vec[i]]->parent==NULL)
                    {
                        (t->children).push_back(ap[vec[i]]);
                        ap[vec[i]]->parent=t;
                        que.push(ap[vec[i]]);
                    }
                }
            }
        }
        return root;
    }

    void order(struct node *root) 
    {
        if(root->isartiPoint)
        {
            cout<<"arti("<<root->artiPoint<<") ";
            for(int i=0;i<(root->children).size();i++)
            {
                order((root->children)[i]);
            }
        }
        else 
        {
            cout<<"bcc(";
            for(int i=0;i<(root->vertices).size();i++)
            {
                cout<<(root->vertices)[i];
            }
            cout<<") ";
            for(int i=0;i<(root->children).size();i++)
            {
                order((root->children)[i]);
            }
        }
    }

    vector<int> findAffectedNodes(vector<Edge> &B)
    {
        queue <int> vertSet;
        vector<int> aff;
        for(int i=0;i<v;i++)
            vertSet.push(i);
        for(int i=0;i<B.size();i++)
        {
            int x=B[i].from;
            int y=B[i].to;
            int siz=vertSet.size();
            for(int j=0;j<siz;j++)
            {
                int w=vertSet.front();
                vertSet.pop();
                if(d[x][w]!=inf && d[y][w]!=inf)
                {
                    if(d[x][w]-d[y][w]>1 || d[x][w]-d[y][w]<-1)
                    {
                        aff.push_back(w);
                        inL[w]=1;
                    }
                    else
                        vertSet.push(w);
                }
                else if (d[x][w]==inf || d[y][w]==inf)
                {
                    aff.push_back(w);
                    inL[w]=1;
                }
                
                else
                    vertSet.push(w);
                
                siz=vertSet.size();
            }
        }
        return aff;
    }

    void finddelta()
    {
        vector<node*> effbcc;
        
        for(int i=0;i<effbcc.size();i++)//hv
        {
            Delta[effbcc[i]]=0;
            for(int j=0;j<(artiInbcc[effbcc[i]]).size();j++)//ai
            {
                delta[bccno[effbcc[i]]][(artiInbcc[effbcc[i]])[j]]=0;
                vector<node*> bccinsubtree;
                subtree((artiInbcc[effbcc[i]])[j] , bccinsubtree);
                for(int k=0;k<(bccinsubtree).size();k++)
                {
                    node* hd =(bccinsubtree)[k];
                    int hdsize=(((bccinsubtree)[k])->vertices).size();
                    int a = (hd->parent)->artiPoint;
                    delta[bccno[effbcc[i]]][(artiInbcc[effbcc[i]])[j]]=delta[bccno[effbcc[i]]][(artiInbcc[effbcc[i]])[j]] + hdsize*(d[(artiInbcc[effbcc[i]])[j]][a]-prevd[(artiInbcc[effbcc[i]])[j]][a]);
                    delta[bccno[effbcc[i]]][(artiInbcc[effbcc[i]])[j]]=delta[bccno[effbcc[i]]][(artiInbcc[effbcc[i]])[j]] + (D[a]-prevD[a]);
                }
                Delta[effbcc[i]] = Delta[effbcc[i]] + delta[ bccno[effbcc[i]] ][(artiInbcc[effbcc[i]])[j]];
            }
        }
    }

    void subtree(int ar,vector<node*> &out)
    {
        node* par = ap[ar];
        if((par->children).size()==0)
        {
            return;
        }
        else 
        {
            for(int i=0;i<(par->children).size();i++)
            {
                node * tmp = (par->children)[i];
                out.push_back(tmp);
                for(int j=0;j<(tmp->children).size();j++)
                {
                    subtree(((tmp->children)[j])->artiPoint, out);
                }
            }
        }

    }

    int noofsub(int ar)
    {
        node* par = ap[ar];
        int out=1;
        if((par->children).size()==0)
        {
            return out;
        }
        else 
        {
            for(int j=0;j<(par->children).size();j++)
            {
                node* tmp=(par->children)[j];
                out++;
                for(int i=0;i<(tmp->children).size();i++)
                {
                    out = out + noofsub(((tmp->children)[i])->artiPoint);
                }
            }
        }
        return out;
    }

    void bottomUpTraversal(node* itr) 
    {
        stack<node*> sta1;
        stack<node*> sta2;
        node* root=BCT();
        sta1.push(root);
        while(!sta1.empty())
        {
            node* tmp=sta1.top();
            sta1.pop();
            sta2.push(tmp);
            for(int i=0;i<(tmp->children).size();i++)
            {
                node * arti = (tmp->children)[i];
                for(int j=0;j<(arti->children).size();j++)
                {
                    node * temp=(arti->children)[j];
                    sta1.push(temp);
                }
            }
        }

        while(!sta2.empty())
        {
            // for(int i=0;i<a.size();i++)
            // {
                node* pb=sta2.top();
                sta2.pop();
                node* gp;
                if(pb->parent != NULL)
                {
                    gp=pb->parent;
                }
                else 
                {
                    gp=NULL;
                }
                for(int j=0;j<(artiInbcc[pb]).size();j++)//v
                {
                    int v=(artiInbcc[pb])[j];
                    if(gp==NULL)
                    {
                        int wt=v-0;
                        
                        // for(int k=0;k<(artiInbcc[pb]).size();k++)
                        // {
                        //     if(j!=k && gp!=ap[artiInbcc[pb][k]])
                        //     {
                        //         int ai= artiInbcc[pb][k];
                        //         delte[bccno[pb]][v]= delte[bccno[pb]][v] + delte[bccno[pb]][ai];
                        //     }
                        // }
                        // delte[bccno[pb]][v] = delte[bccno[pb]][v] + DD[bccno[pb]][v] - D[bccno[pb]][v];
                        // delte[bccno[pb]][v] = delte[bccno[pb]][v] + wt * ( dd[gp][v] - d[gp][v] );

                        // what should we do in last line if gp == NULL ?? ;
                    }
                    else if(artiInbcc[pb][j]!=gp->artiPoint)
                    {
                        int wt=v-noofsub(gp->artiPoint);
                        //int v=artiInbcc[pb][j];
                        for(int k=0;k<(artiInbcc[pb]).size();k++)
                        {
                            if(j!=k && gp!=ap[artiInbcc[pb][k]])
                            {
                                int ai= artiInbcc[pb][k];
                                delta[bccno[pb]][v]= delta[bccno[pb]][v] + delta[bccno[pb]][ai];
                            }
                        }
                        delta[bccno[pb]][v] = delta[bccno[pb]][v] + D[v] - prevD[v];
                        delta[bccno[pb]][v] = delta[bccno[pb]][v] + wt * ( d[gp->artiPoint][v] - prevd[gp->artiPoint][v] );
                    }
                }
            //}
        }
    }

    void topDownTraversal(node* itr) 
    {
        queue<vector<int> > que;
        node* root=BCT();
        que.push(artiInbcc[root]);
        while(!que.empty())
        {
            vector<int> a=que.front();
            que.pop();
            for(int i=0;i<a.size();i++)
            {
                node* pb=(ap[a[i]])->parent;
                node* gp;
                if(pb->parent != NULL)
                {
                    gp=pb->parent;
                }
                else 
                {
                    gp=NULL;
                }
                for(int j=0;j<(artiInbcc[pb]).size();j++)//v
                {
                    int v=artiInbcc[pb][j];
                    if(gp==NULL)
                    {
                        int wt=v-0;
                        
                        // for(int k=0;k<(artiInbcc[pb]).size();k++)
                        // {
                        //     if(j!=k && gp!=ap[artiInbcc[pb][k]])
                        //     {
                        //         int ai= artiInbcc[pb][k];
                        //         delte[bccno[pb]][v]= delte[bccno[pb]][v] + delte[bccno[pb]][ai];
                        //     }
                        // }
                        // delte[bccno[pb]][v] = delte[bccno[pb]][v] + DD[bccno[pb]][v] - D[bccno[pb]][v];
                        // delte[bccno[pb]][v] = delte[bccno[pb]][v] + wt * ( dd[gp][v] - d[gp][v] );

                        // what should we do in last line if gp == NULL ?? ;
                    }
                    else if(artiInbcc[pb][j]!=gp->artiPoint)
                    {
                        int wt=v-noofsub(gp->artiPoint);
                        //int v=artiInbcc[pb][j];
                        for(int k=0;k<(artiInbcc[pb]).size();k++)
                        {
                            if(j!=k && gp!=ap[artiInbcc[pb][k]])
                            {
                                int ai= artiInbcc[pb][k];
                                delta[bccno[pb]][v]= delta[bccno[pb]][v] + delta[bccno[pb]][ai];
                            }
                        }
                        delta[bccno[pb]][v] = delta[bccno[pb]][v] + D[v] - prevD[v];
                        delta[bccno[pb]][v] = delta[bccno[pb]][v] + wt * ( d[gp->artiPoint][v] - prevd[gp->artiPoint][v] );
                    }
                    for(int k=0;k<(ap[v]->children).size();k++)
                    {
                        vector<int> vec=artiInbcc[(ap[v]->children)[k]] ;
                        que.push(vec);
                    }
                }
            }
        }
    }

    vector<double> UpdateBatch (vector<Edge> &B)
    {
        // L = findAffectedNodes(G,B)
        vector<int> L;
        L = findAffectedNodes(B);
        cout<<"L -> ";
        for(auto ll:L)
            cout<<ll<<" ";
        cout<<endl;
        
        // G′=G∪B,
        for(int i=0;i<B.size();i++)
        {
            int x=B[i].from;
            int y=B[i].to;
            Adj[x][y]=1;
            Adj[y][x]=1;
        }

        //printing AdjMatrix
        cout<<endl<<"Adj Matrix :"<<endl;
        for(int i=0;i<v;i++)
        {
            for(int j=0;j<v;j++)
                cout<<Adj[i][j]<<" ";
            cout<<endl;
        }

        // T=BCT(G′)
        findBCC();
        struct node* root = BCT();
        //checking tree
        //cout<<"HI"<<endl;
        order(root);
        //cout<<"HIEE"<<endl;

        //update degree values
        findDegrees();

        // cc′(v)=cc(v)for each v∈G′ : already copied through the default copy constructor

        // R = R3∪R4∪Q in G′
        vector<int> R;
        for(int ii=0;ii<v;ii++)
        {
            if(isR3(ii)||isR4(ii)||isChainNode(ii))
            {
                R.push_back(ii);
                inR[ii]=1;
            }
        }

        cout<<"Entering line 6"<<endl;
        // for each articulation point a ∈ L,  Compute D'[a] through BFS(a, Ha )
        for(int i=0;i<L.size();i++)
        {
            cout<<"i = "<<i<<endl;
            int a=L[i];
            if(articularPoints[a]==0)
                continue;
            double DDD=0;
            int c=0;
            for(auto tt:BCCofatri[a])
            {
                c++;
                cout<<"c = "<<c<<endl;
                DDD+=BFS_Req(a,tt->vertices);
                cout<<"DDD = "<<DDD<<endl;
            }
            cout<<"HELL"<<endl;
            D[a]=DDD/c;
            cout<<"D["<<a<<"] = "<<D[a]<<endl;
        }
        //cout<<"Line 6 done"<<endl;

        int DDD;
        int inS[v];
        vector<int> S_R;

        // for each affected BCC S ∈ G ′
        for(int i=0;i<L.size();i++)
        {
            int u=L[i];
            if(articularPoints[u]==1)
                continue;
            node* S=BCCofVert[u];//bcc.size()
            
            for(int j=0;j<v;j++)
                inS[j]=0;
            for(int j=0;j<S->vertices.size();j++)
            {
                if(inR[S->vertices[j]]==0)
                    S_R.push_back(S->vertices[j]);
                inS[S->vertices[j]]=1;
            }

            // Compute DV′ [v] through BFS(v, S) for all v such that v ∈ S and v ∈ L ∩ R
            if(inS[u]==1 && inR[u]==1 && inL[u]==1)
            {
                DDD=BFS_Req(u,S->vertices);
            }

            // Compute DV′ [v] through BFS(v, S \ R) for all v such that v ∈ S\R and v∈L
            if(inS[u]==1 && inR[u]==0 && inL[u]==1)
            {
                DDD=BFS_Req(u,S_R);
            }
            D[u]=DDD/(BCCofatri[u]).size();
        }

        prevD=D;




        finddelta();
        // bottomUpTraversal(T )// Similar to Algorithm 4 
        bottomUpTraversal(rootBCT);
        // topDownTraversal(T )// Refer Algorithm 4
        topDownTraversal(rootBCT);

        // foreachBCChinT do
        // ∆[ht] = Σi∆h[ai] for all articulation points ai ∈ h
        for(int i=0;i<BCC.size();i++)
        {
            double temp=0;
            vector<int> Vtemp=artiInbcc[BCC[i]];
            for(int j=0;j<Vtemp.size();j++)
            {
                temp=temp+ delta [bccno[BCC[i]]] [Vtemp[j]];
            }
            Delta[BCC[i]]=temp;
        }


        //update cc value
        for(int i=0;i<v;i++)
        {
            if(articularPoints[i]==0)
                CC[i]=(v-1)/D[i]-Delta[BCCofVert[i]];
            else
            {
                double temp=0;
                int c=0;
                for(auto tt:BCCofatri[i])
                {
                    c++;
                    temp+=Delta[BCCofVert[i]];
                }
                temp=temp/c;
                CC[i]=(v-1)/D[i]-temp;
            }

        }

        cout<<endl<<"reached end"<<endl;
        return CC;
    }
};

int main()
{
    int v=6;
    //cin >> v;
    ifstream inp;
    inp.open("input.csv");
    ofstream outfile;
    outfile.open("output.csv");
    outfile<<"graph_no"<<",";
    for(int i=0;i<v;i++)
    {
        outfile<<i<<',';
    }
    outfile<<endl;
    string line;
    getline(inp,line);
    stringstream s(line);
    int n;string wordn;
    getline(s,wordn,',');
    n=stoi(wordn);
    vector<Graph*>all_g(n+1);
    for(int i=0;i<=n;i++)
    {
        all_g[i]=new Graph(v);
    }
    do
    {
        int init_g,final_g;
        stringstream s(line);
        string word;
        getline(s,word,',');
        init_g=stoi(word);
        getline(s,word,',');
        final_g=stoi(word);
        outfile<<final_g<<',';
        (*all_g[final_g])=(*all_g[init_g]);//uses default copy assignment
        vector < Edge > Batch;
        while (getline(s,word,','))
        {
            int temp=stoi(word);
            Edge e1;
            e1.from=temp/v;
            e1.to=temp%v;
            Batch.push_back(e1);
        }
        vector<double>cc=(*all_g[final_g]).UpdateBatch(Batch);
        // cc.push_back(10.5);
        // cc.push_back(235);
        // cc.push_back(10.45);
        // cc.push_back(445.5);

        for(auto c:cc)
        {
            outfile<<c<<',';
        }
        outfile<<endl;
    }
    while(getline(inp,line));
    
}