# my_track
オプションにファイル名や生成するパターン、生成数などを指定することができる。
引数なしで実行するとヘルプが出力される。

マルチスレッドでの処理をできるようにしたい。。。

##出力されるファイル
1. {filename}_tot.dat
	- Time Over Threshold の飛跡データを2 x 1024 x 256 の形で出力される
	- anode, cathode の順番で出力される。
1. {filename}_flush.dat
	- ナマの波形データの形で出力される
	- tot ファイルと同じ形で出力される
1. {filename}_idealvalue.dat
	- 物理的な情報を出力する
	- E [MeV], m [MeV], phi [rad], theta [rad] の順番で出力される
	- 複数の粒子がある場合は上の組が粒子の数だけ続く
1. {filename}_teachervalue.dat
	- range [mm], phi [rad], theta [rad], sca_a_x, sca_a_y, sca_c_x, sca_c_y, end_a_x, end_a_y, end_c_x, end_c_y [pixel] 
	の順に出力される
	- 複数の粒子がある場合は上の首が粒子の数だけ続く
