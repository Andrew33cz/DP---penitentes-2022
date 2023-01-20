/* ------------------------------------
*      bc. Ondrej Cech (xcecho06), 2022
*  ------------------------------------
*/

#include <iostream>
#include <math.h>
#include <vector>

#include "plycpp.h"
#include <filesystem>
#include <array>
#include <map>
#include <algorithm>

#include<random>
#include<chrono>

#include <omp.h>

using namespace std;


/*
 * @brief   struktura pro popis bodu
*/
struct Point3f {
  float x;
  float y;
  float z;
};

/*
 * @brief   hash tabulka, uklada body se stejnymi souradnicemi
*/
multimap<std::string, int> mapped={};

/*
-------POMOCNE FUNKCE--------------
*/

/*
 * @brief   skalarni soucin
 * @param   A vektor
 * @param   B vektor
 * @return  skalarni soucin vektoru A a B
*/
float Dot(Point3f A, Point3f B){
    return A.x*B.x+A.y*B.y+A.z*B.z;
}

/*
 * @brief   delka vektoru
 * @param   vect vektor
 * @return  delka vektoru vect
*/
float length(Point3f vect){
    return sqrt(vect.x*vect.x + vect.y*vect.y + vect.z*vect.z);
}

/*
 * @brief   vypocet delta x
 * @param   src zkoumany bod
 * @param   vis viditelny bod
 * @return  vzdalenost bodu na vodorovne plose
*/
float distXY(Point3f src, Point3f vis){
    return sqrt((src.x-vis.x)*(src.x-vis.x)+(src.y-vis.y)*(src.y-vis.y));
}

/*
 * @brief   vypocet delta h
 * @param   src zkoumany bod
 * @param   vis viditelny bod
 * @return  rozdil vysky bodu
*/
float diffZ(Point3f src, Point3f vis){
    return src.z-vis.z;
}

/*
 * @brief   vypocet h'(x)
 * @param   ray paprsek
 * @param   vis viditelny bod
 * @param   visNorm normala v miste viditelneho bodu
 * @return  derivace v bode x
*/
float derviaceH(Point3f ray, Point3f vis, Point3f visNorm){
    if(ray.x==0.0 && ray.y==0.0){
        //nejde urcit rez, ale nemelo by nastat
        cerr<<"svisly paprsek, nemelo by nastat"<<endl;
        return -1;
    }
    if(visNorm.z==0){
        //kolmy povrch
        //TODO
        return 65535.0;
    }
    float d=-(visNorm.x*vis.x+visNorm.y*vis.y+visNorm.z*vis.z);
    float movx=vis.x+ray.x;
    float movy=vis.y+ray.y;

    float movz=(visNorm.x*movx+visNorm.y*movy+d)/visNorm.z;
    return (vis.z-movz)/sqrt((movx-vis.x)*(movx-vis.x)+(movy-vis.y)*(movy-vis.y));
}

/*
 * @brief   porovnavani vektoru, pokud jsou dva normalizovane vektory stejne, bude vysledek 0
 * @param   A vektor
 * @param   B vektor
 * @return  rozdil A a B
*/
float diff(Point3f a, Point3f b){
    return abs(a.x-b.x+a.y-b.y+a.z-b.z);
}

/*
 * @brief   normalizace vektoru
 * @param   vect vektor
 * @return  vektor stejneho smeru jako vect o delce 1
*/
Point3f normalize(Point3f vect){
    Point3f hlp;
    float l=length(vect);
    if(l==0.0){
        hlp={0.0,0.0,0.0};
    }else{
        hlp.x=vect.x/l;
        hlp.y=vect.y/l;
        hlp.z=vect.z/l;
    }
    return hlp;
}

/*
 * @brief   vektorovy soucin
 * @param   mid bod
 * @param   a bod
 * @param   b bod
 * @return  vektorovy soucin vektoru mid->A a mid->b
*/
Point3f cross(Point3f mid, Point3f a, Point3f b){
    Point3f hlp;
    float a_x= a.x-mid.x;
    float a_y= a.y-mid.y;
    float a_z= a.z-mid.z;

    float b_x= b.x-mid.x;
    float b_y= b.y-mid.y;
    float b_z= b.z-mid.z;
    //vektorovy soucin=vektor kolmy na trojuhelnik
    hlp.x=a_y*b_z-a_z*b_y;
    hlp.y=a_z*b_x-a_x*b_z;
    hlp.z=a_x*b_y-a_y*b_x;
    return hlp;
}

auto print_node = [](const auto &node) {
    cout << "[" << node.first << "] = " << node.second << '\n';
};

auto print_result = [](auto const &pair) {
    cout << (pair.second ? "inserted: " : "assigned: ");
    print_node(*pair.first);
};

/*
-------PRACE S TABULKOU INDEXU--------------
*/

/*
 * @brief   uklada data do hash tabulky podle stejnych souradnic
 * @param   data data souboru PLY
*/
void indexVerts(plycpp::PLYData data){
    auto vertexElement = data["vertex"];
    const float* ptX = vertexElement->properties["x"]->ptr<float>();
	const float* ptY = vertexElement->properties["y"]->ptr<float>();
    const float* ptZ = vertexElement->properties["z"]->ptr<float>();
    #pragma omp parallel for
    for (size_t i = 0; i < vertexElement->size(); ++i)
			{
			    string _hash=to_string(ptX[i])+to_string(ptY[i])+to_string(ptZ[i]);
			    #pragma omp critical
			    {
				mapped.insert({_hash,i});
			    }
			}
}

/*
 * @brief   najde v hash tabulce data se stejnou souradnici
 * @param   _map hash tabulka
 * @param   key klic do tabulky
 * @return  vektor indexu do
*/
vector<int> findSame(multimap<string, int> _map,string key){
        vector<int> help;
        for (auto itr = _map.begin(); itr != _map.end(); itr++)
            if (itr -> first == key)
               help.push_back(itr -> second);
        return help;
}


/*
-------NORMALY--------------
*/

/*
 * @brief   korekce normal v 3d souboru ply pro jeden face
 * @param   data nactena 3d data
 * @param   triangle index face
*/
void setNormals(plycpp::PLYData data, int triangle){
    const auto& vertexIndicesData = data["face"]->properties["vertex_indices"];
			if (vertexIndicesData && vertexIndicesData->isList)
			{
                //indexy bodu v trojuhelniku
				int i1=vertexIndicesData->at<unsigned int>(3*triangle);
				int i2=vertexIndicesData->at<unsigned int>(3*triangle+1);
				int i3=vertexIndicesData->at<unsigned int>(3*triangle+2);

				Point3f pt1={data["vertex"]->properties["x"]->at<float>(i1), data["vertex"]->properties["y"]->at<float>(i1), data["vertex"]->properties["z"]->at<float>(i1)};
                Point3f pt2={data["vertex"]->properties["x"]->at<float>(i2), data["vertex"]->properties["y"]->at<float>(i2), data["vertex"]->properties["z"]->at<float>(i2)};
                Point3f pt3={data["vertex"]->properties["x"]->at<float>(i3), data["vertex"]->properties["y"]->at<float>(i3), data["vertex"]->properties["z"]->at<float>(i3)};
                //vektorovy soucin
                Point3f vc=normalize(cross(pt1, pt2, pt3));
                //skalarni soucin s puvodni normalou? jestli maji podobny smer
                Point3f vc_orig={data["vertex"]->properties["nx"]->at<float>(i1), data["vertex"]->properties["ny"]->at<float>(i1), data["vertex"]->properties["nz"]->at<float>(i1)};
                float dot=Dot(vc, vc_orig);
                if(dot>0.0){
                    //zapsat do vsech tri bodu trojuhelniku
                    data["vertex"]->properties["nx"]->at<float>(i1)=vc.x;
                    data["vertex"]->properties["nx"]->at<float>(i2)=vc.x;
                    data["vertex"]->properties["nx"]->at<float>(i3)=vc.x;

                    data["vertex"]->properties["ny"]->at<float>(i1)=vc.y;
                    data["vertex"]->properties["ny"]->at<float>(i2)=vc.y;
                    data["vertex"]->properties["ny"]->at<float>(i3)=vc.y;

                    data["vertex"]->properties["nz"]->at<float>(i1)=vc.z;
                    data["vertex"]->properties["nz"]->at<float>(i2)=vc.z;
                    data["vertex"]->properties["nz"]->at<float>(i3)=vc.z;
                }else{
                    //zapsat do vsech tri bodu trojuhelniku
                    data["vertex"]->properties["nx"]->at<float>(i1)=-vc.x;
                    data["vertex"]->properties["nx"]->at<float>(i2)=-vc.x;
                    data["vertex"]->properties["nx"]->at<float>(i3)=-vc.x;

                    data["vertex"]->properties["ny"]->at<float>(i1)=-vc.y;
                    data["vertex"]->properties["ny"]->at<float>(i2)=-vc.y;
                    data["vertex"]->properties["ny"]->at<float>(i3)=-vc.y;

                    data["vertex"]->properties["nz"]->at<float>(i1)=-vc.z;
                    data["vertex"]->properties["nz"]->at<float>(i2)=-vc.z;
                    data["vertex"]->properties["nz"]->at<float>(i3)=-vc.z;
                }

			}
			else
				cout << "No valid list of vertex indices." << endl;
}

/*
-------VIDITELNOST BODU--------------
*/

/*
 * @brief   testuje zda se bod nachazi v trojuhelniku
 * @param   tested testovany bod
 * @param   tri1 bod trojuhelniku
 * @param   tri2 bod trojuhelniku
 * @param   tri3 bod trojuhelniku
 * @return  true, pokud se bod nachazi v trojuhelniku
*/
bool inTriangle(Point3f tested, Point3f tri1, Point3f tri2, Point3f tri3){
    Point3f tri=normalize(cross(tri1, tri2, tri3));
    Point3f ax=normalize(cross(tri1, tested, tri3));
    //pokud ax vyslo [0,0,0] lezi na zkoumane hrane a tedy patri do trojuhelniku
    //pokud maji stejnou normalu, lezi na stejne strane hrany jako 3. bod
    if((abs(diff(tri,ax))<0.001)||(ax.x==0.0 && ax.y==0.0 && ax.z==0.0)){
        //zkoumame dalsi hranu obdobne
        ax=normalize(cross(tri1, tri2, tested));
        if((abs(diff(tri,ax))<0.001)||(ax.x==0.0 && ax.y==0.0 && ax.z==0.0)){
            ax=normalize(cross(tri2, tri3, tested));
            tri=normalize(cross(tri2, tri3, tri1));
            //cout<<"res: "<<ax.x<<", "<<ax.y<<", "<<ax.z<<endl;
           // cout<<diff(tri,ax)<<endl;
            if((abs(diff(tri,ax))<0.001)||(ax.x==0.0 && ax.y==0.0 && ax.z==0.0)){
                return true;
            }else{
                return false;
                }
        }else{
            return false;
            }
    }else{
        return false;
        }
}

/*
 * @brief   vypocet pruniku bodu a plochy
 * https://stackoverflow.com/questions/7168484/3d-line-segment-and-plane-intersection
 * @param   contact bod pruniku - chci zjistit
 * @param   ray paprsek
 * @param   rayOrigin bod pocatku paprsku
 * @param   normal normala plochy
 * @param   coord bod na plose
 * @param   param parametr primky (vzdalenost mezi rayOrigin a contact)
 * @return  true, pokud je bod nalezen, contact a param budou obsahovat hodnoty
*/
bool linePlaneIntersection(Point3f* contact, Point3f ray, Point3f rayOrigin, Point3f normal, Point3f coord, float* param) {
    // get d value
    float d = Dot(normal, coord);
    if (Dot(normal, ray) == 0) {
        return false; // No intersection, the line is parallel to the plane
    }
    // Compute the X value for the directed line ray intersecting the plane
    float dff = (d - Dot(normal, rayOrigin)) / Dot(normal, ray);
    if(dff<0.0){return false;}
    // output contact point
    contact->x = rayOrigin.x + normalize(ray).x*dff; //Make sure your ray vector is normalized
    contact->y = rayOrigin.y + normalize(ray).y*dff;
    contact->z = rayOrigin.z + normalize(ray).z*dff;
    *param=dff;
    return true;
}



/*
-------SIMULACE--------------
*/

/*
 * @brief
 * https://math.stackexchange.com/questions/1585975/how-to-generate-random-points-on-a-sphere
 * @param   number pocet paprsku z jednoho bodu
 * @param   point souradnice bodu
 * @param   data nactena 3d data
 * @return  zmena vysky pro bod point
*/
float computeDiffHeight(int number, Point3f point, vector<int> indexes, plycpp::PLYData data){
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<float> distribution(0.0, 0.34);
    float Iota=0;
    int vsblpnt=0;
        #pragma omp parallel for
        for(int i=0; i<number; i++){
                Point3f result;
                Point3f vec;
                vec.x=distribution(generator);
                vec.y=distribution(generator);
                vec.z=distribution(generator);

                float param_m=10000.0;
                Point3f visible;
                Point3f visnorm;
                bool check =false;
                bool found=false;
                    //pokud neni v poloprostoru normaly
                    //if(vec.z<=(normal.x*vec.x+normal.y*vec.y)/normal.z){}
                //ted mam vektor delky 1 do libovolneho smeru polokoule
                vec=normalize(vec);
               // cout<<"koule: "<<vec.x<<" "<<vec.y<<" "<<vec.z<<endl;


                const auto& vertexIndicesData = data["face"]->properties["vertex_indices"];
                //pro kazdy face privraceny alespon k jednomu z bodu vypocitat prunik s rovinnou,
                for(size_t j=0; j<indexes.size();j++){
                    Point3f point_norm;
                    point_norm.x=data["vertex"]->properties["nx"]->at<float>(indexes[j]);
                    point_norm.y=data["vertex"]->properties["ny"]->at<float>(indexes[j]);
                    point_norm.z=data["vertex"]->properties["nz"]->at<float>(indexes[j]);
                    //cout<<"normala vertexu "<<indexes[j]<<": "<<point_norm.x<<" "<<point_norm.y<<" "<<point_norm.z<<endl;
                    //pokud lze vyslat
                    if(Dot(point_norm, vec)>0){
                        check=true;
                        //muze byt vyslan jen jednou
                        break;
                    }//fi muze byt vyslan
                    //pokud ho nejde vyslat z zadneho vrcholu, zahodit
                    if(check==false) continue;
                }//for indexes
                //pres vsechny faces
                for(size_t k=0; k<vertexIndicesData->size()/3; k++){
                    //cout<<k*3<<endl;
                    Point3f face_norm;
                    Point3f face_base;
                    Point3f face_base2;
                    Point3f face_base3;
                        //cout<<"normala facu "<<j<<": ";
                    face_norm.x=data["vertex"]->properties["nx"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));
                    face_norm.y=data["vertex"]->properties["ny"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));
                    face_norm.z=data["vertex"]->properties["nz"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));

                    face_base.x=data["vertex"]->properties["x"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));
                    face_base.y=data["vertex"]->properties["y"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));
                    face_base.z=data["vertex"]->properties["z"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));

                    face_base2.x=data["vertex"]->properties["x"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+1));
                    face_base2.y=data["vertex"]->properties["y"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+1));
                    face_base2.z=data["vertex"]->properties["z"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+1));

                    face_base3.x=data["vertex"]->properties["x"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+2));
                    face_base3.y=data["vertex"]->properties["y"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+2));
                    face_base3.z=data["vertex"]->properties["z"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+2));
                    //cout<<"face norm: "<<face_norm.x<<" "<<face_norm.y<<" "<<face_norm.z<<endl;
                    float p;
                        //contact, ray, rayorigin, normal coord param
                    if(linePlaneIntersection(&result, vec, point, face_norm, face_base, &p)){
                            //kdyz bude param <0, zahodit (paprsek sel dozadu a narazil na odvracenou stranu)
                        if((inTriangle(result, face_base, face_base2, face_base3))){
                        //v rovine vypocitat jestli se nachazi v trojuhelniku,
                            //kdyz ne continue;
                        //kdyz ano vybrat nejblizsi, ostatni zakryva
                                //cout<<p<<endl;
                                if((p<param_m)&&(p>0)){
                                    param_m=p;
                                    visible=result;
                                    visnorm=face_norm;
                                    found=true;
                                }//fi
                            }//fi inTriangle
                        }//fi linePlanet
                }//for faces

                //ve visible bude proste nejblizsi bod. ten bude bud privraceny=viditelny, nebo odvraceny=sli jsme skrz objekt, chyba
                if((Dot(vec,visnorm)<0)&&(found==true)&&(param_m>0.1)){
                    //vypocitat vzorec ablace, sumovat s mezivysledkem
                    float deltax=distXY(point,visible);
                    float deltah=diffZ(point,visible);
                    float derivh=derviaceH(vec,visible, visnorm);

                    vsblpnt++;
                    #pragma omp critical
                    {
                    Iota=Iota+((deltah-derivh*deltax)/(deltax*deltax+deltah*deltah));
                    }
                }//fi muze dopadnout
                //paprsek co nic nenasel zahodit
    }//for numbers
    float I= -0.5e6/7e9*Iota;
    cout<<"visible points: "<<vsblpnt<<endl;
     return I;
}








/*
 * @brief   krok simulace
 * @param   name jmeno souboru typu ply
 * @param   stepLength delka kroku v sekundach
 * @param   steps pocet kroku
 * @param   rays pocet vyslanych paprsku na jeden bod
*/
void simulate(string name, float stepLength, size_t steps, size_t rays){
    plycpp::PLYData data;
    plycpp::PLYData data2;
    plycpp::load(name, data);
    plycpp::load(name, data2);
    for(size_t s=0; s<steps; s++){
        indexVerts(data);
        cout<<"vertices indexed"<<endl;

        auto vertexElement = data["vertex"];
        auto vertexIndicesData = data["face"];
        vector<int> checkedIndexes;
        for(size_t i = 0; i < vertexElement->size(); i++){
            //vypocitej vysku v kroku
            Point3f pnt={data["vertex"]->properties["x"]->at<float>(i),data["vertex"]->properties["y"]->at<float>(i),data["vertex"]->properties["z"]->at<float>(i)};
            //najit vsechny body se stejnymi souradnicemi
            string key=to_string(pnt.x)+to_string(pnt.y)+to_string(pnt.z);
            vector<int> indexes= findSame(mapped, key);
            //pokud uz byl vertex pocitan s jinym, preskocit
            cout<<"checked: "<<checkedIndexes.size()<<" of "<<vertexElement->size()<<" now:"<<indexes.size()<<endl;
            if(!(any_of(checkedIndexes.begin(), checkedIndexes.end(), [&](const int& elem) { return elem == i; }))){
                // cout<<"indexy: ";
                //vypocitam zmenu vysky v bode
               float differ=computeDiffHeight(rays, pnt, indexes, data);
               //a smazat stary
                mapped.erase(key);
                //zapsat vsechny prepocitane indexy
               for(long long unsigned int a=0; a<indexes.size(); a++){
                    //cout<<"index "<<indexes[a]<<","<<endl;
                    data2["vertex"]->properties["z"]->at<float>(indexes[a])=data["vertex"]->properties["z"]->at<float>(indexes[a])-differ*stepLength;
                    //taky je musí zmìnit v hash tabulce, ulozit novy hash
                    cout<<indexes[a]<<", ";
                }
                //cout<<endl;
                checkedIndexes.insert(checkedIndexes.end(), indexes.begin(), indexes.end());
                cout<<" changed by "<<-differ*stepLength<<endl;
            }
        }
        cout<<"fixing normals"<<endl;
        #pragma omp parallel for
        for(size_t j = 0; j <vertexIndicesData->size(); j++){
            setNormals(data2, j);
        }
        string newname= name.substr(0, name.size()-4);
        newname=newname+"-"+to_string(s)+".ply";
        save(newname, data2, plycpp::FileFormat::ASCII);
        plycpp::load(newname, data);
	}
}


/*
 * @brief
 * https://math.stackexchange.com/questions/1585975/how-to-generate-random-points-on-a-sphere
 * @param   number pocet paprsku z jednoho bodu
 * @param   point souradnice bodu
 * @param   data nactena 3d data
*/
void showVisiblePoints(int number, Point3f point, vector<int> indexes, plycpp::PLYData data){
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<float> distribution(0.0, 0.34);

        for(int i=0; i<number; i++){
                Point3f result;
                Point3f vec;
                vec.x=distribution(generator);
                vec.y=distribution(generator);
                vec.z=distribution(generator);

                float param_m=10000.0;
                Point3f visible;
                Point3f visnorm;
                bool found=false;
                    //pokud neni v poloprostoru normaly
                    //if(vec.z<=(normal.x*vec.x+normal.y*vec.y)/normal.z){}
                //ted mam vektor delky 1 do libovolneho smeru polokoule
                vec=normalize(vec);
               // cout<<"koule: "<<vec.x<<" "<<vec.y<<" "<<vec.z<<endl;


                const auto& vertexIndicesData = data["face"]->properties["vertex_indices"];
                //pro kazdy face privraceny alespon k jednomu z bodu vypocitat prunik s rovinnou,
                for(size_t j=0; j<indexes.size();j++){
                    Point3f point_norm;
                    point_norm.x=data["vertex"]->properties["nx"]->at<float>(indexes[j]);
                    point_norm.y=data["vertex"]->properties["ny"]->at<float>(indexes[j]);
                    point_norm.z=data["vertex"]->properties["nz"]->at<float>(indexes[j]);
                    //cout<<"normala vertexu "<<indexes[j]<<": "<<point_norm.x<<" "<<point_norm.y<<" "<<point_norm.z<<endl;
                    //pokud lze vyslat
                    if(Dot(point_norm, vec)>0){
                        break;
                    }//fi muze byt vyslan
                    //muze byt vyslan jen jednou
                }//for indexes
                //pres vsechny faces
                for(size_t k=0; k<vertexIndicesData->size()/3; k++){
                    Point3f face_norm;
                    Point3f face_base;
                    Point3f face_base2;
                    Point3f face_base3;
                        //cout<<"normala facu "<<j<<": ";
                    face_norm.x=data["vertex"]->properties["nx"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));
                    face_norm.y=data["vertex"]->properties["ny"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));
                    face_norm.z=data["vertex"]->properties["nz"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));

                    face_base.x=data["vertex"]->properties["x"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));
                    face_base.y=data["vertex"]->properties["y"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));
                    face_base.z=data["vertex"]->properties["z"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3));

                    face_base2.x=data["vertex"]->properties["x"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+1));
                    face_base2.y=data["vertex"]->properties["y"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+1));
                    face_base2.z=data["vertex"]->properties["z"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+1));

                    face_base3.x=data["vertex"]->properties["x"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+2));
                    face_base3.y=data["vertex"]->properties["y"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+2));
                    face_base3.z=data["vertex"]->properties["z"]->at<float>(data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3+2));
                    //cout<<"face: "<<data["face"]->properties["vertex_indices"]->at<unsigned int>(k*3)<<endl;
                    float p;
                        //contact, ray, rayorigin, normal coord param
                    if(linePlaneIntersection(&result, vec, point, face_norm, face_base, &p)){
                            //kdyz bude param <0, zahodit (paprsek sel dozadu a narazil na odvracenou stranu)
                        if(inTriangle(result, face_base, face_base2, face_base3)){
                        //v rovine vypocitat jestli se nachazi v trojuhelniku,
                            //kdyz ne continue;
                        //kdyz ano vybrat nejblizsi, ostatni zakryva
                                //cout<<p<<endl;
                               // if(p!=0)cout<<result.x<<" "<<result.y<<" "<<result.z<<endl;
                                if((p<param_m)&&(p>0)){
                                    param_m=p;
                                    visible=result;
                                    visnorm=face_norm;
                                    found=true;
                                }//fi
                            }//fi inTriangle
                        }//fi linePlanet
                }//for faces
                //ve visible bude proste nejblizsi bod. ten bude bud privraceny=viditelny, nebo odvraceny=sli jsme skrz objekt, chyba
                if((Dot(vec,visnorm)<0)&&(found==true)){
                   cout<<visible.x<<" "<<visible.y<<" "<<visible.z<<endl;
                }//fi muze dopadnout
                //paprsek co nic nenasel zahodit
    }//for numbers
}

void getNorms(plycpp::PLYData data){
    for(int x=0; x<data["vertex"]->size(); x++){
      cout<<data["vertex"]->properties["x"]->at<float>(x)+data["vertex"]->properties["nx"]->at<float>(x)*0.1<<" "
      <<data["vertex"]->properties["y"]->at<float>(x)+data["vertex"]->properties["ny"]->at<float>(x)*0.1<<" "
      <<data["vertex"]->properties["z"]->at<float>(x)+data["vertex"]->properties["nz"]->at<float>(x)*0.1<<endl;
    }
}

/*
 * @brief   testovaci funkce
 * @param   name jmeno souboru typu ply
 * @param   rays pocet vyslanych paprsku na jeden bod
*/
void points(string name, size_t rays){
    //index bodu
    int i=20;
    plycpp::PLYData data;
    plycpp::load(name, data);
    indexVerts(data);
    auto vertexElement = data["vertex"];
    auto vertexIndicesData = data["face"];
    vector<int> checkedIndexes;
    Point3f pnt={data["vertex"]->properties["x"]->at<float>(i),data["vertex"]->properties["y"]->at<float>(i),data["vertex"]->properties["z"]->at<float>(i)};
    string key=to_string(pnt.x)+to_string(pnt.y)+to_string(pnt.z);
    vector<int> indexes= findSame(mapped, key);
    getNorms(data);
    //showVisiblePoints(rays, pnt, indexes, data);
}
