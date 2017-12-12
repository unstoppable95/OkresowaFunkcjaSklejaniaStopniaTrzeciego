unit Projekt_EAN;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls , Vcl.ToolWin, Vcl.ComCtrls, MyVarType,IntervalArithmetic32and64;


type

  TForm1 = class(TForm)
    Label1: TLabel;
    Label2: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Edit1: TEdit;
    Edit2: TEdit;
    Edit3: TEdit;
    Edit4: TEdit;
    Button1: TButton;
    Button2: TButton;
    Button3: TButton;
    Label6: TLabel;
    Memo1: TMemo;
    Label3: TLabel;
    Label7: TLabel;
    Button6: TButton;
    Label8: TLabel;
    Edit5: TEdit;
    Label9: TLabel;
    Edit6: TEdit;
    Edit7: TEdit;
    Label10: TLabel;
    Label11: TLabel;
    Label12: TLabel;
    Label13: TLabel;
    Edit8: TEdit;
    Edit9: TEdit;
    Edit10: TEdit;
    Edit11: TEdit;
    Button4: TButton;
    Button5: TButton;
    Button7: TButton;
    Label14: TLabel;
    Label15: TLabel;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure Button6Click(Sender: TObject);
    procedure Button7Click(Sender: TObject);
    procedure Button5Click(Sender: TObject);
    procedure Button4Click(Sender: TObject);
 

  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;
   indeks_tab, indeks_tab1 : Integer;
   tab_val_F :  vector;
   tab_val_X :  vector;
    przed_x,przed_f ,przed_x1,przed_f1 : vectorprz;
    przed_xx,przed_xx1, przed_xx2: String;
    przedzialowy_xx :Interval;  // zmienna XX przedzialowa

   var liczba_wezlow,liczba_wezlow1, kod_bledu ,kod_bledu2,kod_bledu3,kod_bledu4: Integer;
   punkt_xx : Extended;
   wynik : Extended;
   wynik_wspl : matrix;
   wynik_wspl_przed ,wynik_wspl_przed1: 	matrixprz;
    pom: String;

    { Private declarations }
implementation

{$R *.dfm}
{$I PSPLVAL.pas}
{$I PSPLCNS.pas}
{$I PSPLVAL_PRZEDZIAL.pas}
{$I PSPLCNS_PRZEDZIAL.pas}

//zatwierdz ite elementy zmiennopozycyjna
procedure TForm1.Button1Click(Sender: TObject);

begin



if indeks_tab <= liczba_wezlow then
   begin
   Val(Edit3.Text, tab_val_X [indeks_tab], kod_bledu);
    //przedzialowa
    str(tab_val_X[indeks_tab], pom) ;
      przed_x[indeks_tab]:=int_read(pom);



   Val(Edit4.Text, tab_val_F [indeks_tab], kod_bledu);
   //przedzialowa
     str(tab_val_F[indeks_tab], pom) ;
      przed_f[indeks_tab]:=int_read(pom);
      Memo1.Lines.Append ('Wczytano dane ! ' ) ;

      if indeks_tab=liczba_wezlow then
      begin
      ShowMessage('Wczytano wszystkie dane !');
      end;

    if indeks_tab<=liczba_wezlow-1 then
    begin
    indeks_tab:=indeks_tab+1;

    end;


   end ;
 end;

//zatwierdz dane zwykla
procedure TForm1.Button3Click(Sender: TObject);
begin
    Val(Edit1.Text, liczba_wezlow, kod_bledu);
     SetLength(tab_val_X, liczba_wezlow+1);
     SetLength(tab_val_F, liczba_wezlow+1);
     //przedzilowa
     SetLength(przed_x,liczba_wezlow+1);
     SetLength(przed_f,liczba_wezlow+1);

     Val(Edit2.Text, punkt_xx, kod_bledu);
       str(punkt_xx,przed_xx);
     SetLength(wynik_wspl,4,liczba_wezlow+1);
       Memo1.Lines.Append ('Wczytano dane ! ' ) ;
end;


 //oblicz przedzialowa
procedure TForm1.Button4Click(Sender: TObject);
var
interval_wynik1                    : Interval;
lewy2, prawy2 ,lewy1 , prawy1: String;
   I ,J: Integer;
begin
 interval_wynik1:=periodsplinevalue_przed(liczba_wezlow1,przed_x1,przed_f1,przedzialowy_xx,kod_bledu3);

      if kod_bledu3<>0 then
      begin
            Memo1.Lines.Clear();
            ShowMessage('Wystapi� b��d !');
          Memo1.Lines.Append('WYST�PI� B��D - NIE MO�NA OBLICZY� WARTO�CI ! ');
          Memo1.Lines.Append('St = ' + IntToStr(kod_bledu3));
      end;

        if kod_bledu3=0 then
        begin
          Memo1.Lines.Clear();
    iends_to_strings(interval_wynik1,lewy1,prawy1);
     Memo1.Lines.Append(' Arytmetyka przedzialowa wynik : ');
     Memo1.Lines.Append('St = ' + IntToStr(kod_bledu3));
      Memo1.Lines.Append ('Wynik wynosi od : ' + lewy1+ ' do : ' + prawy1 );
      end;

        periodsplinecoeffns_przed(liczba_wezlow1,przed_x1,przed_f1,wynik_wspl_przed1,kod_bledu4);
         if kod_bledu4<>0 then
      begin
            Memo1.Lines.Clear();
            ShowMessage('Wystapi� b��d !');
          Memo1.Lines.Append('WYST�PI� B��D - NIE MO�NA OBLICZY� WSPӣCZYNNIK�W ! ');
          Memo1.Lines.Append('St = ' + IntToStr(kod_bledu3));
      end;

        if kod_bledu4=0 then
        begin
             Memo1.Lines.Append('Arytmetyka przedzialowa  wynik  wspolczynniki : ');
           Memo1.Lines.Append('St = ' + IntToStr(kod_bledu4));
            for J := 0 to liczba_wezlow1-1 do
            begin
             for I := 0 to 3 do

            begin
             iends_to_strings(wynik_wspl_przed1[I][J],lewy2,prawy2);
                Memo1.Lines.Append ('a'+'['+ IntToStr(I) + ']' +  '[' + IntToStr(J) + ']' + '=' + ' od ' + lewy2 + ' do ' + prawy2 );
             end;
        end;

        end;
end;

//zatwierdz ite elementy przedzialowa
procedure TForm1.Button5Click(Sender: TObject);
begin
      if indeks_tab1 <= liczba_wezlow1 then
          begin

              przed_x1[indeks_tab1].a:=left_read( Edit8.Text);
              przed_x1[indeks_tab1].b:=right_read( Edit9.Text);
              przed_f1[indeks_tab1].a:=left_read( Edit10.Text);
              przed_f1[indeks_tab1].b:=right_read(Edit11.Text);
              Memo1.Lines.Append ('Wczytano dane ! ' ) ;

              if indeks_tab1=liczba_wezlow1 then
               begin
                ShowMessage('Wczytano wszystkie dane !');
                end;

               if indeks_tab1<=liczba_wezlow1-1 then
                begin
                indeks_tab1:=indeks_tab1+1;
                 end;


          end;

end;

// oblicz przedzialowa z punktu
procedure TForm1.Button6Click(Sender: TObject);
var
interval_wynik                    : Interval;
lewy, prawy : String;
    I ,J: Integer;
begin

  interval_wynik:=periodsplinevalue_przed(liczba_wezlow,przed_x,przed_f,int_read(przed_xx),kod_bledu3);
  //wypisanie warto�ci dla arytmetyki przedzia�owej

  if kod_bledu3<>0 then
      begin
          Memo1.Lines.Clear();
          ShowMessage('Wystapi� b��d !');
           Memo1.Lines.Append('WYST�PI� B��D - NIE MO�NA OBLICZY� WARTO�CI ! ');
          Memo1.Lines.Append('St = ' + IntToStr(kod_bledu3));
      end;

   if kod_bledu3=0 then
     begin
     Memo1.Lines.Clear();
     iends_to_strings(interval_wynik,lewy,prawy);
     Memo1.Lines.Append(' Arytmetyka przedzialowa z punktu wynik : ');
     Memo1.Lines.Append('St = ' + IntToStr(kod_bledu3));
     Memo1.Lines.Append ('Wynik wynosi od: ' + lewy+ ' do : ' + prawy );
      end;

      periodsplinecoeffns_przed(liczba_wezlow,przed_x,przed_f,wynik_wspl_przed,kod_bledu4);

          if kod_bledu4<>0 then
      begin
          Memo1.Lines.Clear();
          ShowMessage('Wystapi� b��d !');
           Memo1.Lines.Append('WYST�PI� B��D - NIE MO�NA OBLICZY� WSPӣCZYNNIK�W! ');
          Memo1.Lines.Append('St = ' + IntToStr(kod_bledu3));
      end;

      if kod_bledu4=0 then
        begin
        Memo1.Lines.Append('Arytmetyka przedzialowa z punktu wynik  wspolczynniki : ');
         Memo1.Lines.Append('St = ' + IntToStr(kod_bledu4));

           for J := 0 to liczba_wezlow-1 do
            begin
             for I := 0 to 3 do

            begin
              iends_to_strings(wynik_wspl_przed[I][J],lewy,prawy);
               Memo1.Lines.Append ('a'+'['+ IntToStr(I) + ']' +  '[' + IntToStr(J) + ']' + '=' + ' od ' + lewy + ' do ' + prawy );
             end;
            end;
        end;


end;


//zatwierdz dane przedzialowa
procedure TForm1.Button7Click(Sender: TObject);
begin
   Val(Edit5.Text, liczba_wezlow1, kod_bledu);
     //przedzilowa
     SetLength(przed_x1,liczba_wezlow1+1);
     SetLength(przed_f1,liczba_wezlow1+1);


    przedzialowy_xx.a:=left_read(Edit6.Text);
    przedzialowy_xx.b:=right_read(Edit7.Text);
      Memo1.Lines.Append ('Wczytano dane ! ' ) ;
end;


//oblicz arytmetyka zwykla
procedure TForm1.Button2Click(Sender: TObject);
var
I ,J: Integer;
begin

      Memo1.Lines.Clear();
     wynik:=periodsplinevalue(liczba_wezlow,tab_val_X,tab_val_F,punkt_xx,kod_bledu2);
      if kod_bledu2<>0 then
      begin
           Memo1.Lines.Clear();
           ShowMessage('Wystapi� b��d !');
          Memo1.Lines.Append('WYST�PI� B��D - NIE MO�NA OBLICZY� WARTO�CI! ');
          Memo1.Lines.Append('St = ' + IntToStr(kod_bledu2));
      end;

      if kod_bledu2=0 then
       begin
         Memo1.Lines.Append('Arytmetyka zmiennopozycyjna wynik : ');
          Memo1.Lines.Append('St = ' + IntToStr(kod_bledu2));
      Memo1.Lines.Append ('Dla podanych danych wartosc to : ' + FloatToStrF(wynik, ffExponent, 15,10 ));

        end;

      periodsplinecoeffns(liczba_wezlow,tab_val_X,tab_val_F, wynik_wspl,kod_bledu2);
     if kod_bledu2<>0 then
      begin
           Memo1.Lines.Clear();
           ShowMessage('Wystapi� b��d !');
          Memo1.Lines.Append('WYST�PI� B��D - NIE MO�NA OBLICZY� WSPӣCZYNNIK�W!');
          Memo1.Lines.Append('St = ' + IntToStr(kod_bledu2));
      end;
      if kod_bledu2=0 then
       begin
        Memo1.Lines.Append('Arytmetyka zmiennopozycyjna wspolczynniki wyniki : ');
       Memo1.Lines.Append('St = ' + IntToStr(kod_bledu2));

       for J := 0 to liczba_wezlow-1 do

             begin
           for I := 0 to 3 do

              begin
               Memo1.Lines.Append ('Wartosc'+'['+ IntToStr(I)+']' + '[' +IntToStr(J)+']' + ' = ' + FloatToStrF(wynik_wspl[I][J], ffExponent, 15,10 ));
             end;
               end;
        end;
end;
end.
