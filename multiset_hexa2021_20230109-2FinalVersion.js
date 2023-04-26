//
//   Cell Emergent-Ishida Model 2022Model 
//   2021.10.18
//   2022. 3.31
//   2023. 3.31      Copy rigth : Takeshi Ishida 
//
var iterate  = 0  ;       // Iteration numeber ;繰り返し数
var cellSize = 8 ;        // Cell spacing (pixel)　;セルの間隔（ピクセル）
var meshNum  =100 ;       // Number of vertical and horizontal cells　;縦横のセル数
var boxNum   =100000;     // Number of polymerized molecules (number of virtual boxes)　;重合分子数（ボックス数）
var molNumMax =100 ;      // Number of molecules polymerized per box　;1ボックス当たりの重合する分子数
var w   = 0.680 ;         // Morphological parameter　w ;形状パラメータｗ 0.65
var mol = 30 ;            // Number of types of molecules　;分子種の数
var x,y ;

var cells   = new Array();   // Cell State　セルの状態
var cells_x = new Array();   // x-coordinate of cell セルのｘ座標
var cells_y = new Array();   // y-coordinate of cell　セルのｙ座標
var cells_boxnum = new Array();   // Box number of polymerized molecules present in each cell　セルに存在する重合分子の番号
var nextCells = new Array(); // For storing the next cell state　次のセル状態の保存
 
var adjacent_cells_x = new Array();   // x sequence number of adjacent cell　隣接セルのｘ配列番号
var adjacent_cells_y = new Array();   // y sequence number of adjacent cell　隣接セルのｙ配列番号
var residualRate = new Array();       // Residual rate by molecular species　各分子の残存率

var box = new Array();                // Polymerized molecule = virtual box 重合分子
var box_num = new Array();            // Number of polymerized molecules of each cell　各セルの重合分子の数
var box_num_pre = new Array();        // Number of polymerized molecules of each cell　(previous time step)　各セルの重合分子の数（1つ前のステップ）
var box_residualRate = new Array();   // Residual rate per polymerized molecule 重合分子の残存率
var box_composition = new Array();    // Molecular configuration within the polymerized molecule 重合分子の構成

var p_entropy = new Array();   // Entropy production エントロピー生成量
var m_entropy = new Array();   // Morphological entropy 形態エントロピー量
var max_m_entropy = new Array();   // Maximum morphological entropy = Configuration entropy when molecules are uniformly distributed 形態エントロピーの最大値（分子が一様に分布しているときのエントロピー）
var total_p_entropy = new Array();   // Time series record of entropy production エントロピー生成量の時系列記録
var total_m_entropy = new Array();   // Time series record of morphological entropy 形態エントロピー量の時系列記録

var mol_c1 = new Array();   //  Number of non-contiguous molecules 6 in polymerized molecule 1 ;重合分子1中で連続していない分子６の数
var mol_c2 = new Array();　　　//  Number of two consecutive molecules 6 in polymerized molecule 1 ;重合分子1中で分子６が2つ連続している数
var mol_c3 = new Array();   //  Number of three consecutive molecules 6 in polymerized molecule 1 ;重合分子1中で分子６が3つ連続している数
var mol_c4 = new Array();   //  Number of four consecutive molecules 6 in polymerized molecule 1 ;重合分子1中で分子６が4つ連続している数
var mol_c5 = new Array();   //  Number of five consecutive molecules 6 in polymerized molecule 1 ;重合分子1中で分子６が5つ連続している数
var mol_c_total = new Array();   // mol_c1 + mol_c2 + mol_c3 + mol_c4 +mol_c5 

var info_product = new Array();    // Production of polymerized molecule 1 per time step ;タイムステップごとの重合分子1の生産量
var info_decomp  = new Array();    // Degradation of polymerized molecule 1 per time step ;タイムステップごとの重合分子1の分解量
var info_total   = new Array();    // Total number of polymerized molecule 1 per time step ;タイムステップごとの重合分子1の総量


var canvas;
var ctx;

var buttonStart;
var buttonRandom;
var buttonReset;
var timer1;
var running = false;
 
window.onload = function()
{
    canvas = document.getElementById('Cell Emergent Model');
    ctx = canvas.getContext('2d');
	// Initialization 初期設定
    initCells();
	// Button settings ボタン設定
    buttonStart  = document.getElementById('buttonStart');
    buttonRandom = document.getElementById('buttonRandom');
    buttonReset  = document.getElementById('buttonReset');

    buttonStart.addEventListener('click', onStart, false);
    buttonRandom.addEventListener('click', randomCells, false);
    buttonReset.addEventListener('click', initCells, false);
    canvas.addEventListener('click', canvasClick, false);
};
 
// Processing of "Start" button  開始ボタン
function onStart(){
    if(running){
        clearInterval(timer1);
        buttonStart.value = "Start";
        running = false;
    } else {
		// Iterative calculation 繰り返し計算
        nextGeneration();
		
        timer1 = setInterval("nextGeneration()", 1);
        buttonStart.value = "Stop";
        running = true;
    }
}
 
// Processing of "Initialization" button or ”Reset” button 初期化、リセットボタン
function initCells(){
	var ni, nj ;
	
    ctx.fillStyle = 'rgb( 180, 180, 180)';   // Background color 背景色
    ctx.fillRect(0,0, canvas.width, canvas.height);
	
	// Preparation of array variables 配列変数の準備
    for(i=0;i<=meshNum;i++){
        cells[i]   = new Array();
        cells_x[i] = new Array();
        cells_y[i] = new Array();
        cells_boxnum[i] = new Array();
        nextCells[i] = new Array();
        adjacent_cells_x[i] = new Array();
        adjacent_cells_y[i] = new Array();
		p_entropy[i] = new Array() ;
        box_num[i]       = new Array();
        box_num_pre[i]   = new Array();
	
	       for(j=0;j<=meshNum;j++){
           cells[i][j]            = new Array();
           cells_boxnum[i][j]     = new Array();
           nextCells[i][j]        = new Array();
           adjacent_cells_x[i][j] = new Array();
           adjacent_cells_y[i][j] = new Array();
		   p_entropy[i][j] = new Array() ;
           box_num[i][j]       = new Array();
           box_num_pre[i][j]   = new Array();

	       for(k=0;k<=2;k++){
              cells_boxnum[i][j][k]  = new Array();
              box_num[i][j][k]       = new Array();
              box_num_pre[i][j][k]   = new Array();
		   }
		}
	}
	
    for(i=1;i<=2;i++){
        box[i]     = new Array();
		box_residualRate[i] = new Array(); 
        box_composition[i]  = new Array();    //
        for(j=0;j<=boxNum;j++){
           box[i][j]     = new Array();
  		   box_residualRate[i][j] = new Array(); 
           box_composition[i][j]  = new Array();    //
	    }
	}
    // Set cell coordinates and initialize the number of cell numerators  セルの座標の設定、セルの分子数の初期化 
    for(i=1;i<=meshNum;i++){
        for(j=1;j<=meshNum;j++){
		 	cells_x[i][j]=(0+((j-1)%2)*cellSize/2+(i-1)*cellSize);
		 	cells_y[i][j]=(0+ (j-1)*cellSize*Math.sqrt(3)/2);
		    p_entropy[i][j] = 0 ;

            for(k=0;k<=mol;k++){
                cells[i][j][k]  = 0;
 			}
        }
    }
	// Initial placement of empty virtual box  空ボックスの初期配置
 	for (i = 1; i <= 2; i++) {
 	for (j = 0; j <= boxNum-1; j++) {
  		box[i][j][0]=0;              // Box state number, 0: empty, 1: polymer state  ボックスの状態番号、０：空、１：高分子状態 
  		box[i][j][1]=(Math.floor(j / meshNum)% meshNum)+1 ;   // x-coordinate of the box ボックスのx座標 
  		box[i][j][2]=(j % meshNum)+1 ;                        // y-coordinate of the box ボックスのy座標 
  		box[i][j][3]=0;              // Box state number, 0: empty, 1: 重合子の亜種 

	}
	}
 	for (i = 1; i <= 2; i++) {
 	    for (j = 0; j <= boxNum-1; j++) {
		 	for (k = 0; k <= 7; k++) {
		 	   box_residualRate[i][j][k] = 0.1;         // Residual rate by orientation in each box  各ボックスのボックスごと、方向別の残存率 
		 	}
	    }
	}  
 	for (i = 1; i <= 2; i++) {
 	    for (j = 0; j <= boxNum-1; j++) {
		 	for (k = 0; k <= molNumMax-1; k++) {
               box_composition[i][j][k]=0;              // Initialization of the composition of each box 各ボックスの組成の初期化
		 	}
		}
	}
	// Relation settings of adjacent cells  隣接セルの設定
    for(i=1;i<=meshNum;i++){
        for(j=1;j<=meshNum;j++){
		 		if ((j-1)%2==1) {
		 		  // 1 :Central 中央
		 		  ni= i;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][0]=ni ;
		 		  adjacent_cells_y[i][j][0]=nj ;

		 		  // 2 :Upper 上
		 		  ni= i;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][1]=ni ;
		 		  adjacent_cells_y[i][j][1]=nj ;

		 		  // 3 :Right 右
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][3]=ni ;
		 		  adjacent_cells_y[i][j][3]=nj ;

		 		  // 4 :Lower 下
		 		  ni= i;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][5]=ni ;
		 		  adjacent_cells_y[i][j][5]=nj ;

		 		  // 5 :Left 左
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][6]=ni ;
		 		  adjacent_cells_y[i][j][6]=nj ;

		 		  // 6 : Upper right 右上
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][2]=ni ;
		 		  adjacent_cells_y[i][j][2]=nj ;

		 		  // 7 : Lower right 右下
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][4]=ni ;
		 		  adjacent_cells_y[i][j][4]=nj ;

		 		  // 8 : Lower left 左下
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;

		 		  // 9 : Upper left 左上
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;
		 		}

		 		if ((j-1)%2==0) {
		 		  // 1 : Central 中央
		 		  ni= i;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][0]=ni ;
		 		  adjacent_cells_y[i][j][0]=nj ;

		 		  // 2 : Upper 上
		 		  ni= i;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][2]=ni ;
		 		  adjacent_cells_y[i][j][2]=nj ;

		 		  // 3 : Right 右
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][3]=ni ;
		 		  adjacent_cells_y[i][j][3]=nj ;

		 		  // 4 : Lower 下
		 		  ni= i;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][4]=ni ;
		 		  adjacent_cells_y[i][j][4]=nj ;

		 		  // 5 : Left 左
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= j;
		 		  adjacent_cells_x[i][j][6]=ni ;
		 		  adjacent_cells_y[i][j][6]=nj ;

		 		  // 6 : Upper right 右上
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;

		 		  // 7 : Lower right 右下
		 		  ni= (i+1)-(Math.floor(i/meshNum))*meshNum;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;

		 		  // 8 : Lower left 左下
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= (j+1)-(Math.floor(j/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][5]=ni ;
		 		  adjacent_cells_y[i][j][5]=nj ;

		 		  // 9 : Upper right 左上
		 		  ni= i+(meshNum-1)-(Math.floor((i-1+meshNum-1)/meshNum))*meshNum;
		 		  nj= j+(meshNum-1)-(Math.floor((j-1+meshNum-1)/meshNum))*meshNum;
		 		  adjacent_cells_x[i][j][1]=ni ;
		 		  adjacent_cells_y[i][j][1]=nj ;
		 		}
        }
    }
	
	// Setting of residual rate of each molecule 各分子の残存率の設定
	residualRate[0] = 0.0;
	residualRate[1] = 0.0;
	residualRate[2] = 0.75; // 0.75
	residualRate[3] = 0.05; // 0.05
	residualRate[4] = residualRate[2];
	residualRate[5] = residualRate[3];
	residualRate[6] = 0.0;
	residualRate[7] = 0.0;
	residualRate[8] = 1.0;
	residualRate[9] = 1.0;
	residualRate[10] = 1.0;
	residualRate[11] = 1.0;
	residualRate[12] = 0.0;
	residualRate[13] = 1.0;
	residualRate[14] = 1.0;
	residualRate[15] = 1.0;
    // Energy resource molecule	
	residualRate[20] = 0.4;
	residualRate[21] = 0.4;  // not used
	// Energy carrier molecule
	residualRate[25] = 0.2;
	residualRate[26] = 0.2;  // not used

    // Graphic Drawing 描画
    redraw();
}
 
// Processing of "Initialization" button  イニシャルボタン；特定のセル設定
function randomCells(){

    //  Initial molecular configuration 分子の初期配置
    for(x=1;x<=meshNum;x++){
        for(y=1;y<=meshNum;y++){
           cells[x][y][0] = 1 ;
           cells[x][y][1] = 1000000 ;
           cells[x][y][6] = 100000*0.5;
           cells[x][y][7] = 100000*0.5;
           cells[x][y][12] = 100000 ;
           cells[x][y][13] = 0;
           cells[x][y][14] = 3;  //3 
           cells[x][y][15] = 9;  //9

           cells[x][y][19] = 0;
           cells[x][y][20] = 100;  // Energy resource 100
           cells[x][y][21] = 0;    // not used
           cells[x][y][25] = 100;  // Energy carrier molecule 100
           cells[x][y][26] = 0;    // not used
           // Only the central part of the computation space
		   if ((x>45)&&(x<58)) { 	
           if ((y>44)&&(y<58)) { 	
           cells[x][y][20] = 1000;  // Energy resource 1000
		   cells[x][y][25] = 500;   // Energy carrier molecule 100
		   }
		   }
           // Representative box number of each lattice cell
           cells_boxnum[x][y][1][0] = 0;
           cells_boxnum[x][y][2][0] = 0;
		   // Number of boxes in each lattice
		   box_num[x][y][1] = 0 ;
		   box_num[x][y][2] = 0 ;
 		   box_num_pre[x][y][1] = 0 ;
		   box_num_pre[x][y][2] = 0 ;
        }
    }
	
	// Initial configuration of molecule 2 and molecule 3　分子２および分子３の初期配置 Case.1
    for(var m=2;m<=3;m++){
    x=50; y=50 ; cells[x][y][m] = 100 ; 
    x=51; y=50 ; cells[x][y][m] = 100 ; 
    x=52; y=50 ; cells[x][y][m] = 100 ; 
    x=53; y=50 ; cells[x][y][m] = 100 ; 
	
    x=50; y=51 ; cells[x][y][m] = 100 ; 
    x=51; y=51 ; cells[x][y][m] = 100 ; 
    x=52; y=51 ; cells[x][y][m] = 100 ; 
    x=53; y=51 ; cells[x][y][m] = 100 ; 
	
    x=49; y=52 ; cells[x][y][m] = 100 ; 
    x=50; y=52 ; cells[x][y][m] = 100 ; 
    x=51; y=52 ; cells[x][y][m] = 100 ; 
    x=52; y=52 ; cells[x][y][m] = 100 ; 
	}

	// Initial configuration of molecules in polymerization molecule 1
//	initial =[7,7,6,7,6,7,6,7,6,7,7,6,7,6,6,7,6,6,7,6,6,7,6,6,7,7,6,6,7,6,6,6,7,6,6,6,7,6,6,6,7,6,6,6,7,6,6,6,6,7,6,6,6,6,7,6,6,6,6,7,6,6,6,6,7,6,6,6,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7];
//	initial =[7,7,7,7,7,7,7,7,7,7,7,6,7,6,6,7,6,6,7,6,6,7,6,6,7,7,6,6,7,6,6,6,7,6,6,6,7,6,6,6,7,6,6,6,7,6,6,6,6,7,6,6,6,6,7,6,6,6,6,7,6,6,6,6,7,6,6,6,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7];
	initial =[6,7,6,7,6,7,6,7,6,7,7,6,7,6,6,7,6,6,7,6,6,7,6,6,7,7,6,7,7,6,6,6,7,6,6,6,7,6,6,6,7,6,6,6,7,6,6,6,6,7,6,6,6,6,7,6,6,7,6,7,6,6,6,6,7,6,7,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7];

	// Initial configuration of virtual box ボックスの初期配置
	for (var i = 0; i <= boxNum-1 ; i++) {
	
        if ((box[1][i][1]>45)&&(box[1][i][1]<58)) { 	
        if ((box[1][i][2]>44)&&(box[1][i][2]<58)) { 	
		// Polymerization reaction of informant (polymerized molecule 1)　情報子（重合分子１）の重合反応
	    build_polymer(11,i,6,7,w) ;  
        
		box[1][i][0]=1 ;
        box[1][i][3]=11 ;
		for (j = 0; j <= molNumMax-1; j++) {
 	       box_composition[1][i][j]=initial[j];
	    }
		
		let NumArray  = calc_mol_continuity2(i,1,6) ;   //
        var Sum_Numarray = (NumArray[0]+NumArray[1]+NumArray[2]+NumArray[3]+NumArray[4])/100;
		console.log("conti_rate=",box[1][i][1],box[1][i][2],
		               Math.floor(NumArray[0]/Sum_Numarray),
					   Math.floor(NumArray[1]/Sum_Numarray),
					   Math.floor(NumArray[2]/Sum_Numarray),
					   Math.floor(NumArray[3]/Sum_Numarray),
					   Math.floor(NumArray[4]/Sum_Numarray)) ; 
		}
	    }
	}
	// Counting the number of polymerized molecules in each cell　各セルの重合分子数のカウント
	//calc_polmerNum3(1) ;
	//calc_polmerNum3(2) ;
	
	// Graphic Drawing　描画
    redraw();
}
 
// Redraw process　全体を再描画
function redraw(){
    for(x=1;x<=meshNum;x++){
    for(y=1;y<=meshNum;y++){
	    // Display empty cells　空セルの表示
        //ctx.strokeStyle='rgb(0, 0,0)';
        ctx.strokeStyle='rgb(255,255,255)';
        ctx.lineWidth= 0.5;
	    ctx.fillStyle = 'rgb( 255,255,255)';
        ctx.beginPath();
        ctx.arc( cells_x[x][y]+cellSize,cells_y[x][y]+cellSize , cellSize/2, 0 , 2*Math.PI);
	    ctx.fill();
	    ctx.closePath();
	    ctx.stroke();
	}
	}
	

    for(x=1;x<=meshNum;x++){
    for(y=1;y<=meshNum;y++){
	
	    // Draw of polymerized molecule 1　重合分子1の表示
        //if (box_num[x][y][1]>0) drawPolymer(x, y, 1);
        if (cells_boxnum[x][y][1][1]>0) drawPolymer(x, y, 11);
        if (cells_boxnum[x][y][1][2]>0) drawPolymer(x, y, 12);

	    // Draw of molecule　各分子の表示
  	    // Draw of molecule 1　分子１の表示
        //drawCell(x, y, 1);
	    // Draw of molecule 2　分子２の表示
        //drawCell(x, y, 2);
	    // Draw of molecule 3 分子３の表示
        //drawCell(x, y, 3);
	    // Draw of molecule 4 分子４の表示
        //drawCell(x, y, 4);
	    // Draw of molecule 5　分子５の表示
        //drawCell(x, y, 5);
	    // Draw of molecule 8 分子８の表示
        //drawCell(x, y, 8);
	    // Draw of molecule 9　分子９の表示
        //drawCell(x, y, 9);
	    // Draw of molecule 10　分子10の表示
        //drawCell(x, y, 10);
	    // Draw of molecule 11　分子11の表示
        //drawCell(x, y, 11);
	    // Draw of molecule 13　分子13の表示
        
		//if (cells[x][y][25]>1000) drawCell(x, y, 25);
        //if (cells[x][y][26]>1000) drawCell(x, y, 26);

	    // Draw of polymerized molecule 2　重合分子2の表示
        if ((box_num[x][y][2]>0)&&(box_num[x][y][2]<=1)) drawPolymer(x, y, 2);
        if ((box_num[x][y][2]>1)&&(box_num[x][y][2]<=4)) drawPolymer(x, y, 3);
        if ((box_num[x][y][2]>4)&&(box_num[x][y][2]<=6)) drawPolymer(x, y, 4);
        if ((box_num[x][y][2]>6))                        drawPolymer(x, y, 5);
		
        // Display of iteration number　繰り返し数の表示
        ctx.fillStyle = 'rgba(225, 225, 225, 0.05)';        // Fill color (translucent)　塗りつぶす色,半透明
        ctx.fillRect(15, meshNum*cellSize*Math.sqrt(3)/2-30, 65, 20);    // Rectangle Drawing　矩形描画
		
		ctx.font = '20pt Arial';
		ctx.fillStyle = 'rgba(0, 0, 0)';
		ctx.fillText(iterate, 20, meshNum*cellSize*Math.sqrt(3)/2-10);
    }
    }
	// Display legend　凡例の表示
/*       if (mm==2) ctx.fillStyle = 'rgb( 255,182,193)';  // Pale pink　薄紅色
       if (mm==3) ctx.fillStyle = 'rgb( 255, 25,147)';    // pink　ピンク色
       if (mm==4) ctx.fillStyle = 'rgb( 255, 25, 25)';    // Yellow 赤色
       if (mm==5) ctx.fillStyle = 'rgb( 225,  0,225)';    // Purple 紫色
*/
	// polymerized molecule 1 重合子１
	/*
	x= 20;
	y= meshNum*cellSize*Math.sqrt(3)/2+30 ;
       ctx.fillStyle = 'rgb( 180,225,180)';  // Light green 薄い緑

	   ctx.beginPath();
       ctx.arc( x,y , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();

	   ctx.font = '12pt Arial';
	   ctx.fillStyle = 'rgba(0, 0, 0)';
	   ctx.fillText("Polymerized molecule 1", x+15, y+5);	

	// polymerized molecule 2 重合子2
	x= 20;
	y= meshNum*cellSize*Math.sqrt(3)/2+50 ;
       ctx.fillStyle = 'rgb( 255,182,193)';  // Light pink 薄紅色
	   
	   ctx.beginPath();
       ctx.arc( x,y , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();

	   ctx.font = '12pt Arial';
	   ctx.fillStyle = 'rgba(0, 0, 0)';
	   ctx.fillText("Polymerized molecule 2 （Number of molecules in the lattice is less than 1)", x+15, y+5);	

	x= 20;
	y= meshNum*cellSize*Math.sqrt(3)/2+70 ;
       ctx.fillStyle = 'rgb( 255, 25,147)';  // Pink ピンク色
	   
	   ctx.beginPath();
       ctx.arc( x,y , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();

	   ctx.font = '12pt Arial';
	   ctx.fillStyle = 'rgba(0, 0, 0)';
	   ctx.fillText("Polymerized molecule 2 （Number of molecules ; greater than 1 and less than 4)", x+15, y+5);	

	x= 20;
	y= meshNum*cellSize*Math.sqrt(3)/2+90 ;
       ctx.fillStyle = 'rgb( 255, 25, 25)';  // Red 赤色
	   
	   ctx.beginPath();
       ctx.arc( x,y , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();

	   ctx.font = '12pt Arial';
	   ctx.fillStyle = 'rgba(0, 0, 0)';
	   ctx.fillText("Polymerized molecule 2 （Number of molecules ; greater than 4 and less than 6)", x+15, y+5);	

	x= 20;
	y= meshNum*cellSize*Math.sqrt(3)/2+110 ;
       ctx.fillStyle = 'rgb( 225,  0,225)';  // Purple 紫色
	   
	   ctx.beginPath();
       ctx.arc( x,y , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();

	   ctx.font = '12pt Arial';
	   ctx.fillStyle = 'rgba(0, 0, 0)';
	   ctx.fillText("Polymerized molecule 2 （Number of molecules ; greater than 6)", x+15, y+5);	
*/
}
 
// Molecular Drawing 分子の描画
function drawCell(xx, yy, mm){
    if (cells[xx][yy][mm]>0) {
	   //　Color set 色の設定
       if (mm==1) ctx.fillStyle = 'rgb( 255,255,255)'; // White 白
       if (mm==2) ctx.fillStyle = 'rgba(255,224,32,0.5)'; // Yellow 黄色
       if (mm==3) ctx.fillStyle = 'rgba(   0,　32,255,0.5)'; // Blue 青
       if (mm==4) ctx.fillStyle = 'rgba( 255,  0,  0,0.5)'; // Red 赤
       if (mm==5) ctx.fillStyle = 'rgb( 160, 32,255)'; // Purple 紫
       if (mm==6) ctx.fillStyle = 'rgb( 255,208,160)'; // Light pink 薄紅色
       if (mm==7) ctx.fillStyle = 'rgb( 160,128, 60)'; // Brown 茶色
       if (mm==8) ctx.fillStyle = 'rgb(  80,208,255)'; // Light blue 水色
       if (mm==9) ctx.fillStyle = 'rgb(   0,192,  0)'; // Green 緑色
       if (mm==10) ctx.fillStyle = 'rgb(  0,255,  0)'; // Green 緑
       if (mm==11) ctx.fillStyle = 'rgb(255, 96,208)'; // Pink ピンク色
       if (mm==12) ctx.fillStyle = 'rgb(255,224, 32)'; // Yellow 黄色
       
	   ctx.fillStyle = 'rgba(255,224,32,0.2)'; // Yellow 黄色
	   
	   // Cell drawing セルの表示
	   ctx.beginPath();
       ctx.arc( cells_x[xx][yy]+cellSize,cells_y[xx][yy]+cellSize , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();
    }	
}

// Drawing of polymerized molecules 重合分子の描画
function drawPolymer(xx, yy, mm){
	   //　　Color set 色の設定
       if (mm==1) ctx.fillStyle = 'rgb( 180,225,180)';  // Light green 薄い緑
       if (mm==2) ctx.fillStyle = 'rgb( 255,182,193)';  // Light pink 薄紅色
       if (mm==3) ctx.fillStyle = 'rgb( 255, 25,147)';  // Pink ピンク色
       if (mm==4) ctx.fillStyle = 'rgb( 255, 25, 25)';  // Red 赤色
       if (mm==5) ctx.fillStyle = 'rgb( 225,  0,225)';  // Purple 紫色

       if (mm==11) ctx.fillStyle = 'rgb( 180,225,180)';  // Light green 薄い緑
       if (mm==12) ctx.fillStyle = 'rgb( 180,180,225)';  // Light blue 薄い青

       // Cell drawing セルの表示
	   ctx.beginPath();
       ctx.arc( cells_x[xx][yy]+cellSize,cells_y[xx][yy]+cellSize , cellSize/2, 0 , 2*Math.PI);
	   ctx.fill();
	   ctx.closePath();
	   ctx.stroke();
   
} 

// Main Process ; Reaction, polymerization and diffusion processes 反応・重合・拡散プロセス
function nextGeneration(){
    var i ;
    var molrate1 ;
	var x,y ;

	iterate = iterate +1 ;
	//console.log("iterate = ", iterate) ;

    // Initialization of entropy generation　エントロピー生成量の初期化
    total_p_entropy[iterate] = 0;

    // Initialization
    mol_c1[iterate] =0 ; mol_c2[iterate] =0 ; mol_c3[iterate] =0 ; mol_c4[iterate] =0 ; mol_c5[iterate] =0 ;
    mol_c_total[iterate] =0 ;
	
	info_product[iterate]=0 ; info_decomp[iterate]=0; info_total[iterate]=0 ;

    // Energy production process ;エネルギーの生産
	for (x = 1; x <= meshNum; x++) {
	for (y = 1; y <= meshNum; y++) {
		
	
	   let NumArray  = calc_mol_continuity(x,y,1,6) ;   // Count the number of consecutive numerator 6
       var Sum_num_array = NumArray[0] + NumArray[1] + NumArray[2] + NumArray[3] + NumArray[4];	  
	   var num_array1=NumArray[1]/Sum_num_array;
	   var num_array2=NumArray[2]/Sum_num_array;
	   var num_array3=NumArray[3]/Sum_num_array;
	   var num_array4=NumArray[4]/Sum_num_array;
	   
       // Calculate the reaction rate of molecule 20 → molecule 25
       var react_rate = (num_array1*0.0 +num_array2*0.0 + num_array3*0.0 + num_array4*4.0-0.60)*1.0*100 ; 
	   if (react_rate > 100) react_rate=100;
	   
	   // Reaction molecule 20 → molecule 25
	   chem_react2(x,y,20,0,25,0,react_rate,0) ; // 

       // Supplementation of energy molecule 20 ;エネルギー分子の補給	   
       if (cells[x][y][20]<=0)  {
        //cells[x][y][20]=cells[x][y][20]+800; // 55  220
        if (iterate<=4000) cells[x][y][20]=cells[x][y][20]+800; // 55  220
        if (iterate>4000) cells[x][y][20]=cells[x][y][20]+200; // 55  220
        if (iterate>6000) cells[x][y][20]=cells[x][y][20]+100; // 55  220
	   }

       if (cells[x][y][21]<=0)  cells[x][y][21]=0 ; // not used

       // Removal of waste molecule 19s	   
       if (cells[x][y][19]>0)  cells[x][y][19]=0 ; 

       // Record the number of consecutive molecure6's per time step 
       if ((box_num[x][y][1]>0)&&(box_num[x][y][2]>0))  {
           mol_c1[iterate] = mol_c1[iterate] + NumArray[0] ;  
           mol_c2[iterate] = mol_c2[iterate] + NumArray[1] ;  
           mol_c3[iterate] = mol_c3[iterate] + NumArray[2] ;  
           mol_c4[iterate] = mol_c4[iterate] + NumArray[3] ;  
           mol_c5[iterate] = mol_c5[iterate] + NumArray[4] ;  
           //if (iterate==10) console.log(box_composition[1][cells_boxnum[x][y][1][0]]) ;
       }
	}
	}
	
    mol_c_total[iterate] = mol_c1[iterate] + mol_c2[iterate] + mol_c3[iterate] + mol_c4[iterate] + mol_c5[iterate] ;

    // RConverting from real numbers to percentages
    mol_c1[iterate] = Math.round(mol_c1[iterate] /mol_c_total[iterate]*100*Math.pow( 10, 1 ))/Math.pow( 10, 1 ) ;	// 小数点以下２桁で四捨五入
    mol_c2[iterate] = Math.round(mol_c2[iterate] /mol_c_total[iterate]*100*Math.pow( 10, 1 ))/Math.pow( 10, 1 ) ;	// 小数点以下２桁で四捨五入
    mol_c3[iterate] = Math.round(mol_c3[iterate] /mol_c_total[iterate]*100*Math.pow( 10, 1 ))/Math.pow( 10, 1 ) ;	
    mol_c4[iterate] = Math.round(mol_c4[iterate] /mol_c_total[iterate]*100*Math.pow( 10, 1 ))/Math.pow( 10, 1 ) ;	
    mol_c5[iterate] = Math.round(mol_c5[iterate] /mol_c_total[iterate]*100*Math.pow( 10, 1 ))/Math.pow( 10, 1 ) ;	
	if ((iterate % 10)==0) console.log("mol_c=",mol_c1[iterate],mol_c2[iterate],mol_c3[iterate],mol_c4[iterate],mol_c5[iterate],mol_c_total[iterate]);
	
	// Diffusion process of energy-related molecules
	diff(20) ;
	//diff(21) ;
	diff(25) ;
	//diff(26) ;

	// Molecular Diffusion　分子の拡散
	diff(1) ;
	diff(2) ;
	diff(3) ;

	// Molecular Reactions　分子の反応  chem_react(x1,x2,y1,y2,r);  x1 + x2 → y1 + y2  Reaction rate r 反応確率　ｒ
	chem_react(2,0,4,0,0,0,5,1) ; // 5
	chem_react(3,0,5,0,0,0,2,1) ; // 2

	// Molecular Diffusion　分子の拡散
	diff(4) ;
    diff(5) ;
	diff(6) ;
	diff(7) ;
//	diff(10) ;
//	diff(11) ;
	diff(12) ;


	// Molecular Reactions 分子の反応  chem_react(x1,x2,y1,y2,r);  x1 + x2 → y1 + y2  Reaction rate r 反応確率　ｒ
	chem_react(4,1,4,8,0,0,100,1) ;
	chem_react(5,1,5,9,0,0,100,1) ;
	chem_react(4,1,4,11,0,1,100,1) ;

 	chem_react(8,0,10,0,0,1,-1,1) ;  // -1；Reacts in the ratio of molecules 6 in the polymerized molecule 重合分子中の分子６の比率で反応する
	chem_react(9,0,10,0,0,1,-1,1) ;  // -1；Reacts in the ratio of molecules 6 in the polymerized molecule 重合分子中の分子６の比率で反応する

	
	chem_react(8,0,1,0,0,0,100,0) ;
	chem_react(9,0,1,0,0,0,100,0) ;

	chem_react(10,11,1,1,0,0,100,0) ;   // Compare numerators 10 and 11, and the one with the higher number remains. 分子10，11のどちらは多いほうが残る

	// Molecular Polymerization 分子の重合
    //cal_polymer(1) ;	

 	for (i = 0; i <= boxNum-1; i++) {
		
		// Polymerization of informant (polymerized molecule 1) 情報子の重合
		if (box[1][i][0]==0) {
		if ((cells[box[1][i][1]][box[1][i][2]][11]>0)) {  // If molecule 11 exists 分子11が存在するなら
		let NumArray  = calc_molrate(1,i,6,7) ;   // Calculation of the ratio of molecules in a polymerized molecule 分子の比率の計算
		molrate1 = NumArray[0];
		pol1 = NumArray[1];
		pol_num = NumArray[2];

		// Set mutation rate; not used 突然変異のレートを設定
		var mutation_rate =0; // not used
		var s = Math.floor( Math.random()*100);	// not used
		// Mutation duplicates informants of subspecies.; not used 突然変異で違う種類の情報子が複製される
        if (s< mutation_rate) {
			if (pol1 ==11) {pol1=12;}  // not used
			else 
			if (pol1 ==12) {pol1=11;}  // not used
		}

		if (molrate1>0) {
            // Molecular Polymerization 分子の重合；pol_numで指定した重合子と同じものをコピーする；一定の突然変異も指定できる
		    if (pol1==11) num =build_polymer_mutate(11,i,6,7,molrate1,pol_num,0.10) ;  // (Polymerizer type, replication destination box number, component molecule 1, component molecule 2, composition ratio, replication source box number, mutation rate)(重合子種類,複製先ボックス番号,構成分子1,構成分子2,構成比率,複製元のボックス番号,突然変異率)
		    if (pol1==12) num =build_polymer_mutate(12,i,6,7,molrate1,pol_num,0.10) ;  // not used
            
			info_product[iterate] = info_product[iterate] +num ;
			// Molecular Polymerization 分子の重合；特定の比率で新たに重合する
		    //if (pol1==11) build_polymer(11,i,6,7,molrate1) ;
		    //if (pol1==12) build_polymer(12,i,6,7,molrate1) ;
			cells[box[1][i][1]][box[1][i][2]][11]=cells[box[1][i][1]][box[1][i][2]][11]-1;
		}  
		}
		}
	}

    // Molecular Reactions
	chem_react(14,11,14,1,0,0,100,1) ;
    chem_react(15,11,15,13,0,0,100,0) ;
	chem_react(13,11,11,11,0,0,100,1) ;

  
 	for (i = 0; i <= boxNum-1; i++) {
		if ((cells[box[2][i][1]][box[2][i][2]][13]>0)) {  //If molecule 13 exists 分子13が存在するなら
			// Formation of molecules 2 and 3 分子２，３の生成  
	        chem_react2(box[2][i][1],box[2][i][2],11,1,11,2,100,1) ;

	        chem_react2(box[2][i][1],box[2][i][2],11,1,11,3,100,1) ;

	        chem_react2(box[2][i][1],box[2][i][2],14,1,14,2,100,1) ;
	        chem_react2(box[2][i][1],box[2][i][2],14,1,14,3,100,1) ;

	        chem_react2(box[2][i][1],box[2][i][2],13,1,13,2,100,1) ;
	        chem_react2(box[2][i][1],box[2][i][2],13,1,13,3,100,1) ;

       // Polymerization of membrane molecules (polymerized molecule 2)　膜分子（重合分子２）の重合
		if (box[2][i][0]==0) { 
          cells[box[2][i][1]][box[2][i][2]][13]=cells[box[2][i][1]][box[2][i][2]][13]-1; 
          cells[box[2][i][1]][box[2][i][2]][1] =cells[box[2][i][1]][box[2][i][2]][1]+1; 
 		  build_polymer(2,i,12,0,1.00) ;  // Molecular Polymerization　分子の重合
		}
		}
	}

	// Counting the number of polymerized molecules in each cell　各セルの重合分子数のカウント
	calc_polmerNum3(1) ;
	calc_polmerNum3(2) ;

	// Diffusion of polymerized molecules　重合分子の拡散
	diff_polymer(1,0.75) ;//0.1 0.7
	diff_polymer(2,0.75) ;//0.75

    cal_polymer(1) ;	

	calc_polmerNum3(1) ;
	calc_polmerNum3(2) ;
	
    // Calculation of entropy in the diffusion process of polymerized molecules 　重合分子の拡散によるエントロピー計算 	
	for (i = 1; i <= meshNum; i++) {
	for (j = 1; j <= meshNum; j++) {
    	for (k = 1; k <= 2; k++) {
    		if (box_num[i][j][k]>=1) { 
			  total_p_entropy[iterate] = total_p_entropy[iterate]   + ((box_num[i][j][k]-box_num_pre[i][j][k])*(box_num[i][j][k]-box_num_pre[i][j][k]))/(box_num[i][j][k]);
			  //total_p_entropy[iterate] = total_p_entropy[iterate]   + ((box_num[i][j][k]-box_num_pre[i][j])*(box_num[i][j][k]-box_num_pre[i][j][k]))/(box_num[i][j][k]);
            }
        }
	}
	}	

	// Molecular Removal　分子の除去
	chem_react(4,0,1,0,0,0,5,0) ; // 5
	chem_react(5,0,1,0,0,0,5,0) ; // 5

	//chem_react(8,0,1,0,0,0,100,0) ;
	//chem_react(9,0,1,0,0,0,100,0) ;

	chem_react(10,0,1,0,0,0,5,0) ; // 5

	// Degradation of polymerized molecules　重合分子の分解
 	for (i = 0; i <= boxNum-1; i++) {
		// Degradation of informant (polymerized molecule 1)　情報子（重合分子１）の分解
		if ((cells_boxnum[box[2][i][1]][box[2][i][2]][2][0]<=0)) { 		
  	       d= decomp_polymer(1,i,0.020) ;  // 分子の分解 0.010
           // Counting the number of polymerized molecule 1 (informant) decompositions ;情報子の分解数のカウント
           info_decomp[iterate] = info_decomp[iterate] + d;
        }
		// Degradation of membrane molecules (polymerized molecules 2)　膜分子（重合分子２）の分解
		if ((cells[box[2][i][1]][box[2][i][2]][13]<=0)) {  // 10 //If molecule 11 exists within a certain range　分子11が一定の範囲で存在するなら
	       decomp_polymer(2,i,0.02);  // Molecular degradation rate 0.05　分子の分解0.02
		}
	}

    // Molecular Reactions
	chem_react(11,0,1,0,0,0,3,0) ; // 3
	chem_react(13,0,1,0,0,0,75,0) ;//75 50 100


    // Calculation of morphological entropy　形態エントロピーの計算
	for (i = 0; i <= mol; i++) {
        m_entropy[i]=0.0;
        max_m_entropy[i]=0.0;
        }
		
    m_entropy[6] =calc_polmerNum2(1) ;
    m_entropy[12]=calc_polmerNum2(2) ;

	for (i = 1; i <= mol; i++) {
        m_entropy[i]=m_entropy[i] + calc_molNum2(i) ;
	    //console.log(i,calc_molNum2(i),m_entropy[i]) ;
        }
	
	for (k = 1; k <= mol; k++) {
		if (m_entropy[k]>0) { 
            max_m_entropy[k]=(m_entropy[k]* Math.log(m_entropy[k]) -m_entropy[k])-(meshNum*meshNum*((m_entropy[k]/(meshNum*meshNum)* Math.log(m_entropy[k]/(meshNum*meshNum)) -m_entropy[k]/(meshNum*meshNum))));
            m_entropy[k]= m_entropy[k]* Math.log(m_entropy[k]) -m_entropy[k];
            }
	    // console.log(k,m_entropy[k]) ;
		}

	for (i = 1; i <= meshNum; i++) {
	for (j = 1; j <= meshNum; j++) {
    	for (k = 1; k <= mol; k++) {
    		if (cells[i][j][k]>=1) { 
    		m_entropy[k]=m_entropy[k] - (cells[i][j][k] * Math.log(cells[i][j][k])-cells[i][j][k]);
            }
        }
	}
	}
	for (k = 1; k <= mol; k++) {
        m_entropy[0]=m_entropy[0] + m_entropy[k];
        max_m_entropy[0]= max_m_entropy[0] + max_m_entropy[k];
    }
		
	total_m_entropy[iterate] = (max_m_entropy[0]- m_entropy[0]);
	//console.log(total_p_entropy[iterate],total_m_entropy[iterate],max_m_entropy[0],m_entropy[0]) ;
	
	// Drawing　描画
    if ((iterate % 10)==0) redraw();
	// Counting of the total number of polymerized molecules1 (informants) ;情報子の総数のカウント	
	info_total[iterate]=calc_polmerNum3(1) ;	
	// Output entropy value　エントロピーデータの出力
	if (iterate==1000) data_download(iterate) ;   // データを出力する繰り返し数を指定する
	if (iterate==2000) data_download(iterate) ;   // データを出力する繰り返し数を指定する
	if (iterate==4000) data_download(iterate) ;   // データを出力する繰り返し数を指定する
	if (iterate==6000) data_download(iterate) ;   // データを出力する繰り返し数を指定する
	if (iterate==8000) data_download(iterate) ;   // データを出力する繰り返し数を指定する

}
 
// Canvas click 
function canvasClick(e){
    var xx = e.clientX - canvas.offsetLeft;
    var yy = e.clientY - canvas.offsetTop;
    var col = Math.floor(xx / cellSize)+1 ;
    var row = Math.floor(yy / (cellSize*(Math.sqrt(3)/2)))+1 ;
    if (cells[col][row][1]==0)  cells[col][row][1]=1; else cells[col][row][1]=0;
	drawCell(col, row,1);
}


// Chemical reaction　化学反応
function chem_react(p1,p2,q1,q2,nmol,npol,r,ene) {
    var pp1, pp2, qq1,qq2 ;
	var rr , rr_ene;
	var flg =1;
	var waste ;
	var ene_rate, ene_stock;

	pp1=0;
	pp2=0;
	qq1=0;
	qq2=0; 
	
	if (p1>0 ) pp1 =1; else pp1=0;
	if (p2>0 ) pp2 =1; else pp2=0;
	if (q1>0 ) qq1 =1; else qq1=0;
	if (q2>0 ) qq2 =1; else qq2=0;

    ene_rate =0.0 ;
	ene_stock=0 ;
    
	var pol_kind = npol ;
	
	if (npol>10) npol=　Math.floor(npol / 10　)  ;	
	
	
	for(var x=1;x<=meshNum;x++){
        for(var y=1;y<=meshNum;y++){
		
		var f=0 ; 
	    if (pol_kind<10) {if (box_num[x][y][npol]>0) f=1;}
	    if (pol_kind>=10) {if (box[npol][cells_boxnum[x][y][npol][0]][3] == pol_kind) f=1;}
			
			if ((nmol==0) || ((nmol>0) && (cells[x][y][nmol]>0))) {
            if ((npol==0) || ((npol>0) && (f==1))) {
			
         	rr_ene = 1.0 ;
			waste = 0;
			var s = Math.floor( Math.random()*100);

            ene_stock = cells[x][y][25] + cells[x][y][26];
            if (ene_stock!=0) ene_rate = cells[x][y][25]/ene_stock ; else ene_rate =0;
            //if (p1==20) console.log(" ene_rate=",cells[x][y][25],cells[x][y][26],ene_stock,ene_rate,x,y) ; 
            //if (Number.isNaN(ene_stock)) console.log("ene_stock=",x,y,ene_stock,cells[x][y][25],cells[x][y][26]);
	        
			var a = cells[x][y][p1];
 	        if (cells[x][y][p1]>0) {
		    

 	        if (p2>0) {
	           if (cells[x][y][p1]<cells[x][y][p2]) {
				    a=cells[x][y][p1]/(qq1+qq2);
			        if ((ene==1) && (cells[x][y][p1]<ene_stock  ))
					{
                      rr_ene = 1.0 ;						
					} else {
					  if (cells[x][y][p1]==0) rr_ene = 0.0 ; else rr_ene = ene_stock  /cells[x][y][p1] ;
			          //if (cells[x][y][p1]==0) rr_ene = 0.0 ;
					}
					
				   } 
			   else {
				    a=cells[x][y][p2]/(qq1+qq2); 
			        if ((ene==1) && (cells[x][y][p2]<ene_stock  ))
					{
                      rr_ene = 1.0 ;						
					} else {
					  if (cells[x][y][p2]==0) rr_ene = 0.0 ; else rr_ene = ene_stock  /cells[x][y][p2] ;
			          //if (cells[x][y][p2]==0) rr_ene = 0.0 ;
					  //if (p1==14) console.log(" rr_ene=",rr_ene,cells[x][y][25],cells[x][y][p2],x,y) ; 
					}
			   }
	        }	
			//if (p1==14) console.log(" rr_ene=",rr_ene,x,y) ; 
 	        if (ene==0) rr_ene = 1.0 ;
	        //console.log(" a=",a) ; 

            if (r==-1) {
			   if (npol==1) rr=calc_molrate2(npol,x,y,6,7); 
			   if (npol==2) rr=1.0; 
			   
			   if (rr>0) {
	           if (p1>0) {cells[x][y][p1] = cells[x][y][p1]-Math.round(a*qq1*rr + a*qq2)*rr_ene ; waste = waste + Math.round(a*qq1*rr + a*qq2)*rr_ene; }
	           if (p1>0) {cells[x][y][p1] = cells[x][y][p1]-Math.round(a*qq1*(1-rr) + a*qq2 )*rr_ene ; waste = waste +Math.round(a*qq1*(1-rr) + a*qq2 )*rr_ene;}

	           if (q1>0) {cells[x][y][q1] = cells[x][y][q1]+Math.round(a*pp1*rr + a*pp2)*rr_ene  ; waste = waste +Math.round(a*pp1*rr + a*pp2)*rr_ene ;}
	           if (q1>0) {cells[x][y][1] = cells[x][y][1]+Math.round(a*pp1*(1-rr) + a*pp2)*rr_ene  ; waste = waste +Math.round(a*pp1*(1-rr) + a*pp2)*rr_ene;}
			   
               if (ene==1) cells[x][y][25] = cells[x][y][25]  -waste *ene_rate;
               if (ene==1) cells[x][y][26] = cells[x][y][26]  -waste *(1.0-ene_rate);
               if (ene==1) cells[x][y][19] = cells[x][y][19]  +waste ;
               if (cells[x][y][25]<0) cells[x][y][25] = 0 ;
               if (cells[x][y][26]<0) cells[x][y][26] = 0 ;
			   
			   // Counting entropy production associated with a reaction　反応に伴うエントロピー生成量のカウント
			   if (q1==1) flg=-1; else flg = 1;
			   total_p_entropy[iterate] = total_p_entropy[iterate] + flg ;

			   }
			} else if (s<=r) {

	           if (p1>0) {cells[x][y][p1] = cells[x][y][p1]-Math.round(a*qq1 + a*qq2)*rr_ene ; waste = waste + Math.round(a*qq1 + a*qq2)*rr_ene;}
	           if (p2>0) {cells[x][y][p2] = cells[x][y][p2]-Math.round(a*qq1 + a*qq2)*rr_ene ; waste = waste + Math.round(a*qq1 + a*qq2)*rr_ene;}
	           if (q1>0) {cells[x][y][q1] = cells[x][y][q1]+Math.round(a*pp1 + a*pp2)*rr_ene ; waste = waste + Math.round(a*pp1 + a*pp2)*rr_ene;}
	           if (q2>0) {cells[x][y][q2] = cells[x][y][q2]+Math.round(a*pp1 + a*pp2)*rr_ene ; waste = waste + Math.round(a*pp1 + a*pp2)*rr_ene;}
			   
               if (ene==1) cells[x][y][25] = cells[x][y][25]- Math.round(waste *ene_rate);
               if (ene==1) cells[x][y][26] = cells[x][y][26]- Math.round(waste *(1.0-ene_rate));
               if (ene==1) cells[x][y][19] = cells[x][y][19]+ Math.round(waste) ;
               if (cells[x][y][25]<0) cells[x][y][25] = 0 ;
               if (cells[x][y][26]<0) cells[x][y][26] = 0 ;
			   
			   // Counting entropy production associated with a reaction　反応に伴うエントロピー生成量のカウント
 			   if (q1==1) flg=-1; else flg = 1;
			   total_p_entropy[iterate] = total_p_entropy[iterate] + flg ;
			   
	        }
			}
			}
			}

		}
	}

}

// Chemical reaction 2 化学反応(セル単位）
function chem_react2(xx,yy,p1,p2,q1,q2,r,ene) {
    var pp1, pp2, qq1,qq2 ;
	var rr , rr_ene;
	var flg =1;
	var waste;
	var ene_rate, ene_stock;
	
	rr_ene = 1.0 ;
	waste =0 ;
	pp1=0;
	pp2=0;
	qq1=0;
	qq2=0; 
	
	if (p1>0 ) pp1 =1 ;
	if (p2>0 ) pp2 =1 ;
	if (q1>0 ) qq1 =1 ;
	if (q2>0 ) qq2 =1 ;
	
    ene_stock = cells[xx][yy][25] + cells[xx][yy][26];
    if (ene_stock!=0) ene_rate = cells[xx][yy][25]/ene_stock ; else ene_rate =0;
	
    var s = Math.floor( Math.random()*100);
        if (s<=r) {
            
	        var a = cells[xx][yy][p1];
 	        if (cells[xx][yy][p1]>0) {
		    
 	        if (p2>0) {
	           if (cells[xx][yy][p1]<cells[xx][yy][p2]) {
				    a=cells[xx][yy][p1]/(qq1+qq2);
			        if ((ene==1) && (cells[xx][yy][p1]<ene_stock))
					{
                      rr_ene = 1.0 ;						
					} else {
					  if (cells[xx][yy][p1]==0) rr_ene = 0.0 ; else rr_ene = ene_stock/cells[xx][yy][p1] ;
			          //if (cells[xx][yy][p1]==0) rr_ene = 0.0 ;
					}
					
				   } 
			   else {
				    a=cells[xx][yy][p2]/(qq1+qq2); 
			        if ((ene==1) && (cells[xx][yy][p2]<ene_stock))
					{
                      rr_ene = 1.0 ;						
					} else {
					  if (cells[xx][yy][p2]==0) rr_ene = 0.0 ; else rr_ene = ene_stock/cells[xx][yy][p2] ;
			          //if (cells[xx][yy][p2]==0) rr_ene = 0.0 ;
					}
			   }
	        }	
 	        if (ene==0) rr_ene = 1.0 ;
			
	        if (p1>0) {cells[xx][yy][p1] = cells[xx][yy][p1]-Math.round(a*qq1 + a*qq2)*rr_ene; waste =waste + Math.round(a*qq1 + a*qq2)*rr_ene;}
	        if (p2>0) {cells[xx][yy][p2] = cells[xx][yy][p2]-Math.round(a*qq1 + a*qq2)*rr_ene; waste =waste + Math.round(a*qq1 + a*qq2)*rr_ene;}
	        if (q1>0) {cells[xx][yy][q1] = cells[xx][yy][q1]+Math.round(a*pp1 + a*pp2)*rr_ene; waste =waste + Math.round(a*pp1 + a*pp2)*rr_ene;}
	        if (q2>0) {cells[xx][yy][q2] = cells[xx][yy][q2]+Math.round(a*pp1 + a*pp2)*rr_ene; waste =waste + Math.round(a*pp1 + a*pp2)*rr_ene;}
	        
            if (ene==1) cells[xx][yy][25] = cells[xx][yy][25]  -waste*ene_rate ;
            if (ene==1) cells[xx][yy][26] = cells[xx][yy][26]  -waste*(1.0-ene_rate) ;
            if (ene==1) cells[xx][yy][19] = cells[xx][yy][19]  +waste;
            if (cells[xx][yy][25]<0) cells[xx][yy][25] = 0 ;
            if (cells[xx][yy][26]<0) cells[xx][yy][26] = 0 ;
			
			// Counting entropy production associated with a reaction　反応に伴うエントロピー生成量のカウント
			if (q1==1) flg=-1; else flg = 1;
			total_p_entropy[iterate] = total_p_entropy[iterate] + flg ;

			}
		}
}

// Molecular Diffusion　分子の拡散
function diff(s) {
	
    for(var x=1;x<=meshNum;x++){
        for(var y=1;y<=meshNum;y++){
	        nextCells[x][y][s] = 0;
		    //if ((s==26)&&(cells[x][y][26]>10000)) console.log(" ene_rate=",cells[x][y][25],cells[x][y][26],x,y) ; 
		}
	}	
	
	for (var i = 1; i <= meshNum; i++) {
	    for (var j = 1; j <= meshNum; j++) {
            t= cells[i][j][s] ;  
            if (t>0) {
			//  Distribution of molecules to adjacent cells　隣接セルへの分配
            for (var p=0 ; p<=6 ; p++) { 
             	cells[i][j][s] =	cells[i][j][s] - Math.floor(t*(1.0-residualRate[s])/7.0);
 			    nextCells[adjacent_cells_x[i][j][p]][adjacent_cells_y[i][j][p]][s] = nextCells[adjacent_cells_x[i][j][p]][adjacent_cells_y[i][j][p]][s] + Math.floor(t*(1.0-residualRate[s])/7.0);
                } 
	  		    // Distribution of remaining molecules to own cell　自身のセルへの残存する分の分配
    		    cells[i][j][s] =cells[i][j][s] - Math.floor(t*residualRate[s]);
    		   	nextCells[i][j][s] = nextCells[i][j][s]+ Math.floor(t*residualRate[s]);
	               // If the divisor has a remainder, it is distributed by a random number.　割り算で余りがでた場合、乱数により割ふりをしてしまう。
     	  	       t= cells[i][j][s];  

    		   	   while (t>0) {  
                    p = Math.floor( Math.random()*7);
  	  		        if (p<=6) { 
  		    	    cells[i][j][s] = cells[i][j][s] - 1;
   				    nextCells[adjacent_cells_x[i][j][p]][adjacent_cells_y[i][j][p]][s] =nextCells[adjacent_cells_x[i][j][p]][adjacent_cells_y[i][j][p]][s] + 1;

      	            } else {
      		    	cells[i][j][s] = cells[i][j][s] - 1;
      		   	    nextCells[i][j][s] = nextCells[i][j][s]+ 1;
                    }
  	  		       t=t-1 ;
              	  }
			}

        }   //  j_loop
	}       //  i_loop
		 
    for(var x=1;x<=meshNum;x++){
        for(var y=1;y<=meshNum;y++){
			
			// Calculation of entropy production associated with molecular diffusion　分子の拡散に伴うエントロピー生成量の計算
			if (cells[x][y][s]>0) {
                total_p_entropy[iterate] = total_p_entropy[iterate]   + (cells[x][y][s]-nextCells[x][y][s])*(cells[x][y][s]-nextCells[i][j][s])/cells[i][j][s];
			}			
			
	        cells[x][y][s] = nextCells[x][y][s];
		}
	}

}

//   Counting of the ratio of molecules in the polymerized molecule (based on cell coordinates)　高分子中の分子の比率のカウント（セル座標ベース）
function calc_molrate2(boxkind,xx,yy,mol1,mol2) {  
	var n1 ,n2 ;
	var rate ;
	var i,j ;
	n1=0 ;
	n2=0 ;
	rate=0;

	// Measures to speed up the calculation of the ratio of polymerized molecule 1 (calculated from one polymerizer without examining all Boxes)　重合分子１の比率の計算を高速化するための措置（全てのBoxを調べずに一つの重合子から計算）
	i =cells_boxnum[xx][yy][boxkind][0] ;
	if (i>0) {
 	for (j = 0; j <= molNumMax-1; j++) {
 	    if (box_composition[boxkind][i][j]==mol1) n1=n1+1 ;
 	    if (box_composition[boxkind][i][j]==mol2) n2=n2+1 ;
 	}
	}
	
	if ((n1 +n2)==0) {rate=0;} else {rate = n1/(n1 +n2) ;}
	//if ((rate>=1.0)&&(rate<1.4)) console.log(n1,n2,rate) ;
	return rate ;
}

//  Count of the ratio of molecules in the polymerized molecule (based on box number)　高分子中の分子の比率のカウント（ボックス番号ベース）
function calc_molrate(boxkind,boxnumber,mol1,mol2) { 
	var n1 ,n2 ;
	var rate=0 ;
	var xx,yy ;
	var i,j ;
    rate=0;
	n1=0 ;
	n2=0 ;
	xx= box[boxkind][boxnumber][1] ;
	yy= box[boxkind][boxnumber][2] ;
	
	// Measures to speed up the calculation of the ratio of polymerized molecule 1 (calculated from one polymerizer without examining all Boxes)　重合分子１の比率の計算を高速化するための措置（全てのBoxを調べずに一つの重合子から計算）
 	i =cells_boxnum[xx][yy][boxkind][0] ;
	if (i>0) {
 	for (j = 0; j <= molNumMax-1; j++) {
 	    if (box_composition[boxkind][i][j]==mol1) n1=n1+1 ;
 	    if (box_composition[boxkind][i][j]==mol2) n2=n2+1 ;
 	}
	}
	
	//rate = n1/(n1 +n2) ;
	if ((n1 +n2)==0) {rate=0;} else {rate = n1/(n1 +n2) ;}
	//if ((rate>=1.00)&&(rate<1.5)) console.log(xx,yy,n1,n2,rate) ;
	return [rate , box[boxkind][i][3], i]
}

// 高分子中の分子の特定の分子の連続度をカウント（座標ベース）
function calc_mol_continuity(xx,yy,boxkind,mol) { 
	var n1 ,n2 ;
	var rate=0 ;
	var xx,yy ;
	var i,j ;
    var premol =0;
	var mol_continu1 =0;
	var mol_continu2 =0;
	var mol_continu3 =0;
	var mol_continu4 =0;
	var mol_continu5 =0;
    var pre_flg =0;

	//xx= box[boxkind][boxnumber][1] ;
	//yy= box[boxkind][boxnumber][2] ;
	
 	i =cells_boxnum[xx][yy][boxkind][0] ;
	if (i>0) {
	if (box_composition[boxkind][i][0]==mol) { premol=mol; pre_flg=1 ;mol_continu1=1;}  
	//premol= box_composition[boxkind][i][0] ; 
 	for (j = 1; j <= molNumMax-1; j++) {
 	    if (box_composition[boxkind][i][j]==mol)
		{
		   mol_continu1 = mol_continu1+1;
		   if (premol==mol) {
			   pre_flg = pre_flg + 1 ;
			   if (pre_flg==2) { mol_continu2 = mol_continu2+1 ;}
			   if (pre_flg==3) { mol_continu3 = mol_continu3+1 ; mol_continu2 = mol_continu2-1 ;}
			   if (pre_flg==4) { mol_continu4 = mol_continu4+1 ; mol_continu3 = mol_continu3-1 ;}
			   if (pre_flg==5) { mol_continu5 = mol_continu5+1 ; mol_continu4 = mol_continu4-1 ; premol=mol;pre_flg =0 ;}
		   } else {
		       premol=mol; pre_flg =1 ;
		   }
		} else {
   	     premol=0; pre_flg =0 ;
		}
		premol = box_composition[boxkind][i][j] ;
 	}
	}
	mol_continu1 = mol_continu1 -(mol_continu2*2+mol_continu3*3+mol_continu4*4+mol_continu5*5);
	return [mol_continu1,mol_continu2, mol_continu3, mol_continu4, mol_continu5]
}

// 高分子中の分子の特定の分子の連続度をカウント（ボックス番号ベース）
function calc_mol_continuity2(i,boxkind,mol) { 
	var n1 ,n2 ;
	var rate=0 ;
	var xx,yy ;
	var i,j ;
    var premol =0;
	var mol_continu1 =0;
	var mol_continu2 =0;
	var mol_continu3 =0;
	var mol_continu4 =0;
	var mol_continu5 =0;
    var pre_flg =0;

	//xx= box[boxkind][boxnumber][1] ;
	//yy= box[boxkind][boxnumber][2] ;
	
// 	i =cells_boxnum[xx][yy][boxkind][0] ;
	if (i>0) {
	if (box_composition[boxkind][i][0]==mol) { premol=mol; pre_flg=1;mol_continu1=1; }  
	//premol= box_composition[boxkind][i][0] ; 
 	for (j = 1; j <= molNumMax-1; j++) {
 	    if (box_composition[boxkind][i][j]==mol)
		{
		   mol_continu1 = mol_continu1+1;
		   if (premol==mol) {
			   pre_flg = pre_flg + 1 ;
			   if (pre_flg==2) { mol_continu2 = mol_continu2+1 ;}
			   if (pre_flg==3) { mol_continu3 = mol_continu3+1 ; mol_continu2 = mol_continu2-1 ;}
			   if (pre_flg==4) { mol_continu4 = mol_continu4+1 ; mol_continu3 = mol_continu3-1 ;}
			   if (pre_flg==5) { mol_continu5 = mol_continu5+1 ; mol_continu4 = mol_continu4-1 ; premol=mol;pre_flg =0 ;}
		} else {
		   premol=mol; pre_flg =1 ;
		}
		} else {
		premol=0; pre_flg =0 ;
		}
		premol = box_composition[boxkind][i][j] ;
 	}
	}
	mol_continu1 = mol_continu1 -(mol_continu2*2+mol_continu3*3+mol_continu4*4+mol_continu5*5);
	
	return [mol_continu1, mol_continu2, mol_continu3, mol_continu4, mol_continu5]
}

// Molecular Polymerization　分子の重合
function build_polymer (boxkind,boxnumber,mol1,mol2,rate) {
	var i,j,k ;
	var p ;
	var n1 ,n2 ;
	var molrate ;
	var rr , rr_ene;
	var flg =1;
	var waste;
	var ene_rate, ene_stock;
    
	n1=0 ;
	n2=0 ;	
    var pol_kind = boxkind ;
	
	if (boxkind>10) boxkind=　Math.floor(boxkind / 10　)  ;
    

 	if (box[boxkind][boxnumber][0]==0)  { //If the box is empty　ボックスが空である

 	if (cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol1]>molNumMax) {  // If there is a material molecule　材料となる分子がある
 	if ((mol2==0)||((mol2>0)&&(cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol2]>molNumMax))) {  // If there is a material molecule　材料となる分子がある

	// エネルギー消費
	ene_rate=1.0;
	ene_stock = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25] + cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][26];
    if (ene_stock!=0) ene_rate = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25]/ene_stock ;

	if (((boxkind==2)&&((molNumMax*0.1)<=ene_stock))||((boxkind==1)&&((molNumMax*0.1)<=ene_stock))) {
	//if (1==1) {
	//console.log(boxkind)
	
	box[boxkind][boxnumber][0]=1 ;
	box[boxkind][boxnumber][3]=pol_kind ;

    
	// 分子の重合
  	    for (i = 0; i <= molNumMax-1; i++) {
            p = Math.floor( Math.random()*100)+1; 	
			if ((rate*100)<=(p)) {
	           box_composition[boxkind][boxnumber][i]=mol2;                    // Description of composition in box　各ボックスの組成の記述 		
		       cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol2] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol2]-1 ;
               n2=n2+1 ;
		    } else {
	           box_composition[boxkind][boxnumber][i]=mol1;                    // Description of composition in box　各ボックスの組成の記述 		
		       cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol1] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol1]-1 ;
               n1=n1+1 ;
		    }
	    }
        //if (boxkind==2) {
		cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25]  -molNumMax*ene_rate *0.01;
        cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][26] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][26]  -molNumMax*(1.0-ene_rate)*0.01 ;
        cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][19] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][19]  +molNumMax*0.01;
        if (cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25]<0) cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25] = 0 ;
        if (cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][26]<0) cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][26] = 0 ;
		//}
		
	// Calculation of entropy production associated with molecular polymerization　分子の重合に伴うエントロピー生成量の計算
	total_p_entropy[iterate] = total_p_entropy[iterate] + molNumMax ;
    }
	}
	}
	}
	molrate = n1/(n1 +n2) ;
	if ((n1 +n2)==0) molrate=0;
	if ((molrate>0)&&(molrate<0.6)) console.log(boxkind,boxnumber,n1,n2,molrate,rate) ;

}

// Molecular Polymerization　分子の複製
function build_polymer_mutate(boxkind,boxnumber,mol1,mol2,rate,pol_num,mutate_rate) {
	var i ;
	var p ;
	var n1 ,n2 ;
	var molrate ;
	var rr , rr_ene;
	var flg =1;
	var waste;
	var ene_rate, ene_stock;
	var num_product ;
    
	n1=0 ;
	n2=0 ;	
	num_product =0;
    var pol_kind = boxkind ;
	
	if (boxkind>10) boxkind=　Math.floor(boxkind / 10　)  ;
    

 	if (box[boxkind][boxnumber][0]==0)  { //If the box is empty　ボックスが空である

 	if (cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol1]>molNumMax) {  // If there is a material molecule　材料となる分子がある
 	if ((mol2==0)||((mol2>0)&&(cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol2]>molNumMax))) {  // If there is a material molecule　材料となる分子がある

	// エネルギー消費
	ene_rate=1.0;
	ene_stock = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25] + cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][26];
    if (ene_stock!=0) ene_rate = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25]/ene_stock ;

	if (((boxkind==2)&&((molNumMax*0.1)<=ene_stock))||((boxkind==1)&&((molNumMax*0.1)<=ene_stock))) {
	//if (1==1) {
	//console.log(boxkind)
	
	box[boxkind][boxnumber][0]=1 ;
	box[boxkind][boxnumber][3]=pol_kind ;

    
	// 分子の重合
  	    for (i = 0; i <= molNumMax-1; i++) {

            //if (p>(mutate_rate*100)) {
			 box_composition[boxkind][boxnumber][i]=box_composition[boxkind][pol_num][i] ;
			//} else {
			//   if (box_composition[boxkind][pol_num][i]==mol1)  box_composition[boxkind][boxnumber][i]=mol2 ;  
			//   if (box_composition[boxkind][pol_num][i]==mol2)  box_composition[boxkind][boxnumber][i]=mol1 ;  
			//}
    	    if (box_composition[boxkind][boxnumber][i]==mol2) {
		       cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol2] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol2]-1 ;
               n2=n2+1 ;
		    }
			if (box_composition[boxkind][boxnumber][i]==mol1) {                    // Description of composition in box　各ボックスの組成の記述 		
		       cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol1] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][mol1]-1 ;
               n1=n1+1 ;
		    }
	    }

        // 突然変異
        
		var p = Math.floor( Math.random()*100); 	
        if (p<=(mutate_rate*100)) {
            var pp = Math.floor( Math.random()*100); 	
			m = box_composition[boxkind][boxnumber][pp];
            i = pp+1; flg=1 ;
			
			while (flg==1) {
               if (m!=box_composition[boxkind][boxnumber][i]){
			     if (box_composition[boxkind][boxnumber][i]==mol1) {
					box_composition[boxkind][boxnumber][i] =mol2 ;
					box_composition[boxkind][boxnumber][pp]=mol1 ;
				    flg=0 ;
				 } 					
			     if (box_composition[boxkind][boxnumber][i]==mol2) {
					box_composition[boxkind][boxnumber][i] =mol1 ;
					box_composition[boxkind][boxnumber][pp]=mol2 ;
				    flg=0 ;
				 }					
			  //flg=0 ;
			  }
              i=i+1; if(i>99) i=0;  
			}
		}
		
        // 重合に伴うエネルギーの消費
		cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25]  -molNumMax*ene_rate *0.01;
        cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][26] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][26]  -molNumMax*(1.0-ene_rate)*0.01 ;
        cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][19] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][19]  +molNumMax*0.01;
        if (cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25]<0) cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][25] = 0 ;
        if (cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][26]<0) cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][26] = 0 ;
		
		
	// Calculation of entropy production associated with molecular polymerization　分子の重合に伴うエントロピー生成量の計算
	total_p_entropy[iterate] = total_p_entropy[iterate] + molNumMax ;
    num_product = num_product +1 ;
	}
	}
	}
	}
	molrate = n1/(n1 +n2) ;
	if ((n1 +n2)==0) molrate=0;
	//console.log(boxkind,boxnumber,n1,n2,molrate,rate);
	if ((molrate>0)&&(molrate<0.4)) console.log(boxkind,boxnumber,n1,n2,molrate,rate) ;

    return num_product ;
}

// Degradation of polymerized molecules　重合分子の分解
function decomp_polymer (boxkind,boxnumber,decomp_rate) {
	
	var molnum ;
    var pol_kind = boxkind ;
	var num_decom ;
	
	num_decomp =0 ;
	if (boxkind>10) boxkind=　Math.floor(boxkind / 10　)  ;
	
 	if (box[boxkind][boxnumber][0]==1)  { //If the box is not empty　ボックスが空でない

    var p = Math.floor( Math.random()*100)+1; 	
	
	if (p<=(decomp_rate*100)) {
	    box[boxkind][boxnumber][0]=0 ;
	    box[boxkind][boxnumber][3]=0 ;
	
  	    for (var i = 0; i <= molNumMax-1; i++) {
			molnum = box_composition[boxkind][boxnumber][i];
			//if (boxkind==1) console.log("composition=",molnum) ;

	        box_composition[boxkind][boxnumber][i]=0;      // Set the composition of the box to 0 (empty)　各ボックスの組成の記述 		
		    cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][molnum] = cells[box[boxkind][boxnumber][1]][box[boxkind][boxnumber][2]][molnum]+1 ;
	    }
    	// Calculation of entropy production associated with the decomposition of polymerized molecules　重合分子の分解に伴うエントロピー生成量の計算
		total_p_entropy[iterate] = total_p_entropy[iterate] + molNumMax ;
    
	    num_decomp = num_decomp + 1 ;
	}
	}
    return num_decomp ;
}

// Diffusion of polymerized molecules　重合分子の拡散
function diff_polymer(mol,rate) {
	var s ;

	// Measures to speed up the calculation of the ratio of polymerized molecules1　重合分子１の比率の計算を高速化するための措置
/*	if (mol==1) {	
    for(x=1;x<=meshNum;x++){
        for(y=1;y<=meshNum;y++){
           cells_boxnum[x][y][mol][0] = 0;
        }
    }	
	}
*/
    for (var i = 0; i <= boxNum-1; i++) {
       	//if(box[mol][i][0]==1) {
        
		    var p = Math.floor( Math.random()*100)+1; 	
		    var u = Math.floor( Math.random()*7); 	
	  		
			//s = calc_polmerNum(mol,adjacent_cells_x[box[mol][i][1]][box[mol][i][2]][u],adjacent_cells_y[box[mol][i][1]][box[mol][i][2]][u]) ;
			s = box_num[adjacent_cells_x[box[mol][i][1]][box[mol][i][2]][u]][adjacent_cells_y[box[mol][i][1]][box[mol][i][2]][u]][mol];
	  		box_residualRate[mol][i][u] = rate + (1-rate)/0.5 *(1/(1+Math.exp(-1*s))-0.5);

	  		if ((box_residualRate[mol][i][u]*100)<=p) {

	  		    box[mol][i][1]= adjacent_cells_x[box[mol][i][1]][box[mol][i][2]][u] ;
		 	    box[mol][i][2]= adjacent_cells_y[box[mol][i][1]][box[mol][i][2]][u] ;
	  	        
			
			}
	  	//}
	
	// Measures to speed up the calculation of the ratio of polymerized molecules1　重合分子１の比率の計算を高速化するための措置
/*	if (mol==1) {	
	if(box[mol][i][0]==1) cells_boxnum[box[mol][i][1]][box[mol][i][2]][mol][0] = i;
	}
*/
    }  // i_loop
}

// Counting the number of polymerized molecules 高分子の数の算
function cal_polymer(mol) {

	if (mol==1) {	
    
	for(x=1;x<=meshNum;x++){
        for(y=1;y<=meshNum;y++){
           cells_boxnum[x][y][mol][0] = 0;
           cells_boxnum[x][y][mol][1] = 0;
           cells_boxnum[x][y][mol][2] = 0;
        }
    }	
    var s = Math.floor( Math.random()*boxNum);
	
	for (var i = 0; i <= boxNum-1; i++) {
      j= i + s ;
      if (j>(boxNum-1)) j=j-(boxNum-1) ; 

	if(box[mol][j][0]==1) {
		cells_boxnum[box[mol][j][1]][box[mol][j][2]][mol][0] = j;
		if (box[mol][j][3]==11) cells_boxnum[box[mol][j][1]][box[mol][j][2]][mol][1] =cells_boxnum[box[mol][j][1]][box[mol][j][2]][mol][1] + 1;
		if (box[mol][j][3]==12) cells_boxnum[box[mol][j][1]][box[mol][j][2]][mol][2] =cells_boxnum[box[mol][j][1]][box[mol][j][2]][mol][2] + 1;
	
	    }
    }  // 
	
	}
}



// Counting of polymers (per cell)　高分子数のカウント（セル単位）
function calc_polmerNum(boxkind,xx,yy) {  
	var n ;
	var number ;
	n=0 ;
	
 	for (var i = 0; i <= boxNum-1; i++) {
	 	if ((box[boxkind][i][1]==xx)&&(box[boxkind][i][2]==yy)) {
           if (box[boxkind][i][0]==1)  n=n+1 ;
		}
	}
	return n ;
}

// Counting of virtual box (per polymerized molecule)　高分子数のカウント（重合分子単位）
function calc_polmerNum2(boxkind) {  
	var n ;
	var number ;
	n=0 ;
	
 	for (var i = 0; i <= boxNum-1; i++) {
           if (box[boxkind][i][0]==1)  n=n+1 ;
	}
	return n ;
}

// Counting of macromolecules (polymerized molecular units) 高分子数のカウント（重合分子セル単位）
function calc_polmerNum3(boxkind) {  
	var n ;
	var number ;
	n=0 ;
	
	for (var i = 1; i <= meshNum; i++) {
	for (var j = 1; j <= meshNum; j++) {	
 	    box_num_pre[i][j][boxkind] = box_num[i][j][boxkind] ;
        box_num[i][j][boxkind] =0 ;
	}
	}
	
	for (var i = 0; i <= boxNum-1; i++) {
           if (box[boxkind][i][0]==1)  {
			 box_num[box[boxkind][i][1]][box[boxkind][i][2]][boxkind] = box_num[box[boxkind][i][1]][box[boxkind][i][2]][boxkind] +1 ;    
			 n=n+1 ;
	       }
	}
	return n ;
}

// Counting of all molecules　全分子のカウント
function calc_molNum() {  
    var num_particle =0;
		for (var i = 1; i <= meshNum; i++) {
		for (var j = 1; j <= meshNum; j++) {
        for (var s = 1; s <= mol  ; s++) {
       		  num_particle = num_particle +cells[i][j][s];    // Counting of molecules　分子数のカウント
        }
		}   //  j_loop
		}   //  i_loop

	  	for (i = 1; i <= 2; i++) {
	  	for (j = 0; j <= boxNum-1; j++) {
	  	 	if (box[i][j][0]==1){
	 	  		for (var k = 0; k <= molNumMax-1; k++) {
	 	  		  if (box_composition[i][j][k]>0) num_particle = num_particle +1 ;
	 		  	}
	 	 	}
	 	}
		}
    return num_particle ;
}
		
// Counting the number of molecules by molecule type　分子種別の分子数のカウント
function calc_molNum2(molkind) {  
    var num_particle =0;
		for (var i = 1; i <= meshNum; i++) {
		for (var j = 1; j <= meshNum; j++) {
        //if (Number.isNaN(cells[i][j][molkind])) console.log(i,j,molkind,cells[i][j][molkind]) ;
        //if (Number.isNaN(cells[i][j][molkind])) cells[i][j][molkind]=0;
		if (cells[i][j][molkind]>0) {			
       		  num_particle = num_particle +cells[i][j][molkind];  // Counting of molecules　　分子数のカウント
		}
		}   //  j_loop
		}   //  i_loop

    return num_particle ;
}

// CSV output of entropy data　エントロピーデータのCSV出力
function data_download(data_len) {

    var str = "";      // Create empty string　空の文字列を作成

    for(var i = 1; i<=data_len; i++){
        str += i+","+total_p_entropy[i]+","+total_m_entropy[i]+","+mol_c1[i]+","+mol_c2[i]+","+mol_c3[i]+","+mol_c4[i]+","+mol_c5[i]+","+info_product[i]+","+info_decomp[i]+","+info_total[i]+"\n"; // Generate output data　出力データを作成
    }

    var blob =new Blob([str],{type:"text/csv"}); //Set the above string(str) in the array　配列に上記の文字列(str)を設定
    var link =document.createElement('a');
    link.href = URL.createObjectURL(blob); 
    link.download ="tempdate.csv";
    link.click();
	
}
		