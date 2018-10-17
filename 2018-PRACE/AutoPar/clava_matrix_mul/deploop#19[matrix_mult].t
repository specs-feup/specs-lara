!------------------------------------------------- arrays
!------ A_1 -> C
integer A_1(0:9999,0:9999)
!------------------------------------------------- loop indices
integer i,j,l
!------------------------------------------------- variables
integer petit_tmp,N,K,M
!------------------------------------------------- body code
for i  =  0  to  N  do
	for j  =  0  to  K  do
		for l  =  0  to  M  do
		petit_tmp = A_1(i,j)
		A_1(i,j) = petit_tmp
		endfor
	endfor
endfor

!---------------------------------------------------------------------------------------------------
!									 Petit Output Dependency 
!-
!-	anti    12: A_1(i,j)        -->  13: A_1(i,j)        (0,0,+)         [ Mo]	 canbeignored = true
!-	flow    13: A_1(i,j)        -->  12: A_1(i,j)        (0,0,+)         [ Mo]	 canbeignored = true
!-	output  13: A_1(i,j)        -->  13: A_1(i,j)        (0,0,+)         [ Mo]	 canbeignored = true
!---------------------------------------------------------------------------------------------------
