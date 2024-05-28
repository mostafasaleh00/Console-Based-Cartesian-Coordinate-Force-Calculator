#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <istream>



using namespace std;

double Fx, Fy, Fz, Fxtotal = 0, Fytotal = 0, Fztotal = 0, Fprime;
double magnitude, magnitudePV, angle, resultant, theta, thetax, thetay, thetaz, unitV[3];
double Frx, Fry, Frz, adjacent, opposite, hypotenuse;
double alpha, peta, omega;
int choice,choice1, axis, type, missing, missing1, positionV[3], coordinates[2][3], choicemissing, angelorratio,negativeplan;
char axis2, axis3, axis1;
double missingFX, missingFY, missingFZ;
int Fcount , Fdirection;
ofstream file("explanation.txt");


void negative_convert(){
    cin>>axis1;
    switch (axis1){
        case 'x':
            Fx=-1*Fx;
            break;
        case 'y':
            Fy=-1*Fy;
            break;
        case 'z':
            Fz=-1*Fz;
            break;
    }
}

void force_direction() {
cout<<"Press '1' if there is ONLY one force in the negative direction :\n"
      "Press '2' if there is TWO force in the negative direction :\n"
      "Press '3' if all the forces are in the negative direction :\n"
      "Press '4' if all forces are in positive direction :";
cin>>Fdirection;

if (Fdirection==1){
    cout<<"Enter the axis (in lowercase):";
    negative_convert();
}
else if (Fdirection==2){
    for (int i = 0; i < 2; ++i) {
        cout<<"Enter axis "<<i+1<<" (in lower case) :";
        negative_convert();
    }
}
else if (Fdirection==3){
    Fx=-1*Fx;
    Fy=-1*Fy;
    Fz=-1*Fz;
}

}

void coorTcart() {
    cout << "enter the magnitude of the force :";
    cin >> magnitude;
    cout << "enter the 2 points :";

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            cin >> coordinates[i][j];
        }
    }
    for (int k = 0; k < 3; k++) {
        positionV[k] = coordinates[1][k] - coordinates[0][k];
    }
    magnitudePV = sqrt(pow(positionV[0], 2) + pow(positionV[1], 2) + pow(positionV[2], 2));
    for (int i = 0; i < 3; i++) {
        unitV[i] = (positionV[i] / magnitudePV) * magnitude;
        if (i == 0)
            Fx = unitV[0];
        else if (i == 1)
            Fy = unitV[1];
        else if (i == 2)
            Fz = unitV[2];
    }
    force_direction();
    file <<"\nsince F"<<Fcount+1<<" is not in cartesian form we going to convert it from\n"
                                "the coordinate system to the cartesian form :";
    file <<"\n1-we going are to calculate the position vector r by \n"
           "subtracting the start point from the end point : \n"
           "r = ("<<positionV[0]<<","<<positionV[1]<<","<<positionV[2]<<")"<<endl;
    file <<"2-calculate the magnitude of r :\n"
           "sqrt( ("<<positionV[0]<<")^2 + ("<<positionV[1]<<")^2 +("<<positionV[2]<<")^2 )"<<endl;
    file <<"3-find the unit vector by dividing position vector r and its magnitude : "<<endl<<
    "U =  ( ("<<unitV[0]<<") i + ("<<unitV[1]<<") j + ("<<unitV[2]<<") k )"<<endl;
    file <<"4-multiply F" <<Fcount+1<<" by Unit vector U to get all three components of the force :\n";
    file <<"Fx = "<<Fx<<endl<<"Fy = "<<Fy<<endl<<"Fz = "<<Fz<<endl;
}

void doubleprojection() {
    cout << "Enter magnitude :\n";
    cin >> magnitude;
    cout << "Enter the first angel\n(the angel between Fprime and the force) :\n";
    cin >> alpha;
    cout << "Enter the second angel and the two axes\n"
            "(the angle between the Fprime and axis) :";
    cin >> peta >> axis2>>axis1;
    cout << "Enter the third angel and the axis\n"
            "separated by a space :";
    cin >> omega >> axis3;
    Fprime = magnitude * cos(alpha * M_PI/180);

    if (axis2 == 'x' && axis3 == 'z') {
        Fx = Fprime * cos(peta * M_PI / 180);
        Fy = Fprime * sin(peta * M_PI / 180);
        Fz = magnitude * cos(omega * M_PI / 180);
    }
    else if (axis2 == 'y' && axis3 == 'z') {
        Fx = Fprime * sin(peta * M_PI / 180);
        Fy = Fprime * cos(peta * M_PI / 180);
        Fz = magnitude * cos(omega * M_PI / 180);
    }
    else if (axis2 == 'z' && axis3 == 'x') {
        Fz = Fprime * cos(peta * M_PI / 180);
        Fy = Fprime * sin(peta * M_PI / 180);
        Fx = magnitude * cos(omega * M_PI / 180);
    }
    else if (axis2 == 'y' && axis3 == 'x') {
        Fy = Fprime * cos(peta * M_PI / 180);
        Fz = Fprime * sin(peta * M_PI / 180);
        Fx = magnitude * cos(omega * M_PI / 180);
    }
    else if (axis2 == 'x' && axis3 == 'y') {
        Fx = Fprime * cos(peta * M_PI / 180);
        Fz = Fprime * sin(peta * M_PI / 180);
        Fy = magnitude * cos(omega * M_PI / 180);
    }
    else if (axis2 == 'z' && axis3 == 'y') {
        Fz = Fprime * cos(peta * M_PI / 180);
        Fx = Fprime * sin(peta * M_PI / 180);
        Fy = magnitude * cos(omega * M_PI / 180);
    }
    force_direction();
    file <<"since that F"<<Fcount+1<<" is not in cartesian form we are going to convert it \n"
           "using the double projection technique :";
    file <<"\nF' = F"<<Fcount+1<<"cos("<<alpha<<") = "<<Fprime<<endl
    <<"F"<<axis2<<" = "<<Fprime<<" cos("<<peta<<") = "<<Fprime * cos(peta * M_PI / 180)<<endl
    <<"F"<<axis1<<" = "<<Fprime<<" sin("<<peta<<") = "<<Fprime * sin(peta * M_PI / 180)<<endl
    <<"F"<<axis3<<" = "<<Fcount+1<<" cos("<<omega<<") = "<<magnitude * cos(omega * M_PI / 180)<<endl;
}


void planforce() {
    cout << "Enter the magnitude :";
    cin >> magnitude;
    cout << "\nEnter the axes the force acts upon\n"
            "(first axis should be the axis the angle is with) :";
    cin >> axis2 >> axis3;
    cout << "Enter '1' for angle \n"
            "Enter '2' for ratio triangle :";
    cin >> angelorratio;

    file <<"\nsince that the F"<<Fcount+1<< " lies on plan we are going to convert it using this equation :";
    if (angelorratio == 1) {
        cout << "\nEnter the angle :";
        cin >> angle;
        Fx = magnitude * cos(angle * M_PI / 180);
        Fy = magnitude * sin(angle * M_PI / 180);
        file <<"\nF"<<axis2<<" = F"<<Fcount+1<<"cos("<<angle<<") = "<<Fx<<endl<<
        "F"<<axis3<<" = F"<<Fcount+1<< "sin("<<angle<<") = "<<Fy<<endl;
    }
    else if (angelorratio == 2) {
        cout << "\nEnter adjacent , opposite and hypotenuse separated by a space :";
        cin >> adjacent >> opposite >> hypotenuse;
        Fx = magnitude * (adjacent / hypotenuse);
        Fy = magnitude * (opposite / hypotenuse);
        file <<"\nF"<<axis2<<" = F"<<Fcount+1<<"* adjacent / hypotenuse = "<<Fx<<endl<<
             "F"<<axis3<<" = F"<<Fcount+1<< "* opposite / hypotenuse = "<<Fy<<endl;
    }
    if (axis2 == 'x' && axis3 == 'z') {
        Fz = Fy;
        Fy = 0;
    }
    else if (axis2 == 'y' && axis3 == 'x') {
        Fz = Fy;
        Fy = Fx;
        Fx = Fz;
        Fz = 0;
    }
    else if (axis2 == 'y' && axis3 == 'z') {
        Fz = Fy;
        Fy = Fx;
        Fx = 0;
    }
    else if (axis2 == 'z' && axis3 == 'x') {
        Fz = Fx;
        Fx = Fy;
        Fy = 0;
    }
    else if (axis2 == 'z' && axis3 == 'y') {
        Fz = Fx;
        Fx = 0;
    }
    force_direction();
}

void dirTcart() {
    cout << "enter '1' if you have a missing angel\n"
            "'2' if you don't\n";
    cin >> choicemissing;
    cout << "Enter The Magnitude Of The Force :";
    cin >> magnitude;

    if (choicemissing == 1) {
        cout << "Enter the axis the missing angle is with :";
        cin >> axis2;
        file <<"since that F"<<Fcount+1 <<"has a missing angle we going to calculate it \n"
               "using the following equation :\n"
               "cos(alpha)^2+cos(peta)^2+cos(gama)^2 = 1";
        if (axis2 == 'x') {
            cout << "Enter peta then gama separated by a space :";
            cin >> peta >> omega;
            alpha = acos(sqrt(1 - pow(cos(omega * M_PI / 180), 2) - pow(cos(peta * M_PI / 180), 2))) * 180 / M_PI;
            file <<"\nalpha = cos^-1(sqrt(1-cos(peta)^2-cos(gama)^2) = "<<alpha<<endl;
        }
        else if (axis2 == 'y') {
            cout << "Enter alpha then gama separated by a space :";
            cin >> alpha >> omega;
            peta = acos(sqrt(1 - pow(cos(omega * M_PI / 180), 2) - pow(cos(alpha * M_PI / 180), 2))) * 180 / M_PI;
            file <<"\npeta = cos^-1(sqrt(1-cos(alpha)^2-cos(gama)^2) = "<<peta<<endl;
        }
        else if (axis2 == 'z') {
            cout << "Enter alpha then peta separated by a space :";
            cin >> alpha >> peta;
            omega = acos(sqrt(1 - pow(cos(alpha * M_PI / 180), 2) - pow(cos(peta * M_PI / 180), 2))) * 180 / M_PI;
            file <<"\ngama = cos^-1(sqrt(1-cos(peta)^2-cos(alpha)^2) = "<<omega<<endl;
        }
    }
    else if (choicemissing == 2) {
        cout << "Enter The Direction Angel With Space Between Each other :";
        cin >> alpha >> peta >> omega;
    }
    Fx = magnitude * cos(alpha * M_PI / 180);
    Fy = magnitude * cos(peta * M_PI / 180);
    Fz = magnitude * cos(omega * M_PI / 180);
    force_direction();
    file <<"\nwe are to convert the force to cartesian form by using those equations :\n"
           "Fx = F"<<Fcount+1<<"cos(alpha) = "<<Fx<<endl <<
           "Fy = F"<<Fcount+1<<"cos(peta) = "<<Fy<<endl<<
           "Fz = F"<<Fcount+1<<"cos(gama) = "<<Fz<<endl;
}

void threeDtype() {
    cout << "\nChoose What Form Is Your Force Written In\n";
    cout << "Enter '1' For Direction Angels\n";
    cout << "Enter '2' if the force lies on a plan\n";
    cout << "Enter '3' For Double Projection\n";
    cout << "Enter '4' For cartesian Coordinates\n";
    cout << "Enter '5' For cartesian form \n";
    cin >> type;

    if (type == 1) { dirTcart(); }
    else if (type == 2) { planforce(); }
    else if (type == 3) { doubleprojection(); }
    else if (type == 4) { coorTcart(); }
    else if (type == 5) {
        cout << "Enter Fx,Fy and Fz Separated By A Space : ";
        cin >> Fx >> Fy >> Fz;

        file <<"\nsince F"<<Fcount+1<<" is already in cartesian form we don't have to convert it : \n";
    }
    Fxtotal += Fx;
    Fytotal += Fy;
    Fztotal += Fz;
}

void threeD() {
    cout << "Enter '1' To Find The Resultant \n";
    cout << "Enter '2' To Find The Missing Force: \n";

    cin >> missing;
    if (missing == 1) {
        cout << "Choose How Many Forces You Want To Find Its Resultant :";
        cin >> choice;

        for (int i = 0; i < choice; i++) {
            Fcount=i;
            threeDtype();
        }

        resultant = sqrt(pow(Fxtotal, 2) + pow(Fytotal, 2) + pow(Fztotal, 2));
        thetax = acos(Fxtotal / resultant);
        thetay = acos(Fytotal / resultant);
        thetaz = acos(Fztotal / resultant);
        cout << "Resultant = " << resultant << endl;
        cout << "thetax = " << thetax * 180 / M_PI << "\nthetay = " << thetay * 180 / M_PI << "\nthetaz = " << thetaz * 180 / M_PI;

        file <<"\nto calculate the resultant we are going to sum all X , Y and Z components :"
               "\nFr = {("<<Fxtotal<<") i + ("<<Fytotal<<") j +("<<Fztotal<<") k }\n";
        file << "\nthen use this equation to calculate the magnitude of the resultant :"
                "\nFr = sqrt ( ("<<Fxtotal<<")^2 + ("<<Fytotal<<")^2 + ("<<Fztotal<<")^2 ) = "<<resultant<<endl;
        file <<"\nto calculate the direction we are going to use these equations :"
               "\ntheta X = cos^-1("<<Fxtotal<<"/"<<resultant<<") = "<<thetax*180/M_PI<<
               "\ntheta Y = cos^-1("<<Fytotal<<"/"<<resultant<<") = "<<thetay*180/M_PI<<
               "\ntheta Z = cos^-1("<<Fztotal<<"/"<<resultant<<") = "<<thetaz*180/M_PI<<endl;

    }

    else if (missing == 2) {
        cout << "enter the resultant force (Frx,Fry,Frz):";
        cin >> Frx >> Fry >> Frz;
        cout << "how many forces you have :";
        cin >> choice;
        for (int i = 0; i < choice - 1; i++) {
            Fcount=i;
            threeDtype();
        }

        missingFX = Frx - (Fxtotal);
        missingFY = Fry - (Fytotal);
        missingFZ = Frz - (Fztotal);
        magnitude = sqrt(pow(missingFX ,2) + pow(missingFY, 2) + pow(missingFZ, 2));
        cout << "the missing force is : \n{" << missingFX << "i + " <<"("<< missingFY <<")"<< " j +" << "("<<missingFZ<<")k }";
        cout<<"\nmagnitude of the resultant = "<<magnitude;
        alpha = acos(missingFX/ magnitude);
        peta = acos(missingFY / magnitude);
        omega = acos(missingFZ / magnitude);
        cout<<"\nthetax = "<<alpha*180/M_PI<<"\nthetay = "<<peta*180/M_PI<<"\nthetaz = "<<omega*180/M_PI;

        file <<"\nto calculate the missing force we going to sum all X , Y and Z components of the known forces : "
        <<"\nFx = "<<missingFX<<" / Fy = "<<missingFY<<" / Fz = "<<missingFZ<<endl;

        file <<"\nsubtract the sum X , Y and Z components of the known forces from Fr : "
               "\nF"<<Fcount+2<<"x = Frx - Fx = "<<missingFX<<endl
               <<"\nF"<<Fcount+2<<"y = Fry - Fy = "<<missingFY<<endl
               <<"\nF"<<Fcount+2<<"z = Frz - Fz = "<<missingFZ<<endl;

        file << "\nthen use this equation to calculate the magnitude of the missing force :"
                "\nF"<<Fcount+2<< " = sqrt ( ("<<missingFX<<")^2 + ("<<missingFY<<")^2 + ("<<missingFZ<<")^2 ) = "<<magnitude<<endl;

        file <<"\nto calculate the direction we are going to use these equations :"
               "\ntheta X = cos^-1("<<missingFX<<"/"<<magnitude<<") = "<<alpha*180/M_PI<<
             "\ntheta Y = cos^-1("<<missingFY<<"/"<<magnitude<<") = "<<peta*180/M_PI<<
             "\ntheta Z = cos^-1("<<missingFZ<<"/"<<magnitude<<") = "<<omega*180/M_PI<<endl;
    }
}



void fxfycalculator() {
    cout << "Choose What Form Is Your Force Written In \n"
            "Enter '1' For Polar Form \n"
            "Enter '2' For Cartesian Form : \n";
    cin >> type;
    if (type == 1) {
        cout << "Enter Magnitude :";
        cin >> magnitude;
        cout << "if the force with angle enter '1' \nif its with ratio triangle enter '2':";
        cin >> angelorratio;
        file <<"since that F"<<Fcount+1<<" is in polar form we are going \n"
              "to convert it to cartesian form by using this equation:\n";
        if (angelorratio == 1) {

            cout << "For Angle With X axis PRESS 1\n";
            cout << "For Angle With Y axis PRESS 2\n";
            cin >> axis;
            cout << "Enter The angle : \n";
            cin >> angle;
            if (axis == 1){
                angle = angle;
                file <<"Fx = "<<magnitude<<"cos("<<angle<<")"<<"= "<<magnitude * cos(angle * M_PI / 180)<<endl;
                file <<"Fy = "<<magnitude<<"sin("<<angle<<")"<<"= "<<magnitude * sin(angle * M_PI / 180)<<endl;
                }
            else if (axis == 2){
                angle = 90 - angle;
                file <<"Fx = "<<magnitude<<"sin("<<angle<<")"<<"= "<<magnitude * cos(angle * M_PI / 180)<<endl;
                file <<"Fy = "<<magnitude<<"cos("<<angle<<")"<<"= "<<magnitude * sin(angle * M_PI / 180)<<endl;

            }
            Fx = magnitude * cos(angle * M_PI / 180);
            Fy = magnitude * sin(angle * M_PI / 180);
            force_direction();
        }

        else if (angelorratio == 2) {
            cout << "For triangle With X axis PRESS 1\n";
            cout << "For triangle With Y axis PRESS 2\n";
            cin >> axis;
            if (axis == 1) {
                cout << "Enter adjacent , opposite and hypotenuse separated by a space :";
                cin >> adjacent >> opposite >> hypotenuse;

                file<<"Fx = "<<magnitude<<"*"<<adjacent<<"/"<<hypotenuse<<"= "<<magnitude * (adjacent / hypotenuse)
                <<"\nFy = "<<magnitude<<"*"<<opposite<<"/"<<hypotenuse<<"= "<<magnitude * (opposite/ hypotenuse)<<endl;

                Fx = magnitude * (adjacent / hypotenuse);
                Fy = magnitude * (opposite / hypotenuse);
            }
            else if (axis == 2) {
                cout << "Enter adjacent , opposite and hypotenuse separated by a space :";
                cin >> adjacent >> opposite >> hypotenuse;

                file<<"Fx = "<<magnitude<<"*"<<opposite<<"/"<<hypotenuse<<"= "<<magnitude * (opposite / hypotenuse)
                    <<"\nFy = "<<magnitude<<"*"<<adjacent<<"/"<<hypotenuse<<"= "<<magnitude * (adjacent / hypotenuse)<<endl;

                Fy = magnitude * (adjacent / hypotenuse);
                Fx = magnitude * (opposite / hypotenuse);
                force_direction();
            }
        }
    }
    if (type == 2) {

        file<<"since F"<<Fcount+1<<" is already in cartesian form we don't have to convert it"<<endl;

        cout << "Enter Fx and Fy Separated By A Space : ";
        cin >> Fx >> Fy;
    }

}

void twoD() {

    cout << "Enter '1' For Finding Resultant and Theta:\n";
    cout << "Enter '2' For Finding Missing Force\n";
    cin >> missing;

    // finding the resultant and theta
    if (missing == 1) {
        cout << "Choose How Many Forces You Want To Find Its Resultant :";
        cin >> choice;

        for (int i = 0; i < choice; i++) {
            Fcount=i;
            fxfycalculator();
            Fxtotal += Fx;
            Fytotal += Fy;
        }
        resultant = sqrt(pow(Fxtotal, 2) + pow(Fytotal, 2));
        theta = atan(Fytotal / Fxtotal);
        file<<"\nnow we will sum all X components and all Y components :";
        file <<"\nFrx = "<<Fxtotal<<"\n"<<"Fry = "<<Fytotal<<'\n';
        file <<"to calculate the resultant and the direction :\n"
               "resultant = sqrt(Frx^2+Fry^2) = "<<"sqrt("<<Fxtotal<<"^2 + "<<Fytotal<<"^2) = "<<resultant<<endl
               <<"theta = tan^-1(Fry/Frx) = "<<"tan^-1("<<Fytotal<<"/"<<Fxtotal<<") = "<<theta* 180 / M_PI;

        cout << "\nResultant = " << resultant << "\ntheta = " << theta * 180 / M_PI<<endl;
    }

        // finding the force1

    else if (missing == 2) {

        //knowing the resultant
        cout << "Enter '1' If The Resultant in Polar Form\n";
        cout << "Enter '2' If It is in Cartesian Form:";
        cin >> resultant;

        if (resultant == 1) {
            cout << "Enter The Magnitude Of Resultant :";
            cin >> magnitude;
            cout << "Enter The Theta Of The Resultant:\n(theta should be with the x axis)";
            cin >> angle;

            Frx = magnitude * cos(angle * M_PI / 180);
            Fry = magnitude * sin(angle * M_PI / 180);

            file <<"since the resultant is in polar form we are going to convert it \n"
                   "to cartesian form using this equation :\n";
            file <<"Frx = F*cos("<<angle<<") = "<<Frx<<endl<<"Fry = F*sin("<<angle<<") = "<<Fry<<endl;
        }
        else if (resultant == 2) {
            cout << "Enter Frx and Fry Separated By A Space :";
            cin >> Frx >> Fry;
            file <<"since resultant is already in cartesian form we don't have to convert it ";
        }

        cout << "Choose How Many Force You Have :";
        cin >> choice;

        for (int i = 0; i < choice - 1; i++) {
            fxfycalculator();
            missingFX = Frx - (Fxtotal);
            missingFY = Fry - (Fytotal);

        }
        file <<"now we could calculate the missing force using this equation :"<<endl<<
        "the missing Fx = Frx-Fxtotal(the sum of all x components) = "<<missingFX<<endl
        <<"the missing Fy = Fry-Fytotal(the sum of all y components) = "<<missingFY;
        cout << "The Missing Force Is :" << endl << "{ Fx = " << missingFX << "  Fy = " << missingFY << " }";
    }

}




int main() {
    int choice1;
    cout << "==============================================================================\n";
    cout << "====================================Hello!====================================\n";
    cout << "===========Welcome To Our Program For Solving Cartesian Coordinates===========\n";
    cout << "==============================================================================\n";
    cout << "==============================================================================\n";

    cout << "\n\n For Begining PRESS Any Number U LIKE <3\n";
    int x;
    cin >> x;
    cout << "Choose The Cartesian Form ";
    cout << "For 2D Enter '2'" << " " << "For 3D Enter '3'" << endl;
    cin >> choice1;
    if (choice1 == 2) {
        twoD();
    }
    else if (choice1 == 3) {
        threeD();
    }
    file.close();
    cout<<"\nEnter '1' for the steps \n"
          "Enter '2' to exit ";
    cin>>choice;
    if(choice==1) {
        ifstream file("explanation.txt");
        cout << "\nSteps :\n";
        cout << file.rdbuf();
        file.close();
    }
    return 0;
}