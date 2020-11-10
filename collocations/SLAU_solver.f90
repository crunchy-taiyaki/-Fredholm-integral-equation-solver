module SLAU_solver
    use precision
         implicit none
         contains

         subroutine gauss(M,X)
         implicit none
         real(mp),dimension(1:) :: X(:)
         real(mp),dimension(1:size(X)+1,1:size(X)) :: M(:,:)
         integer(4) :: n,i,j,k
         real(mp) :: eps
         character(10) :: method='gauss'
         eps=1.0_mp**(-mp)
         n = size(X, dim=1)         
         do k=1,n-1
            if (abs(M(k,k)) < eps) write(*,*) 'Ведущий элемент близок к нулю, номер строки:', k 
            forall (j=k:n+1) M(j,k)=M(j,k)/M(k,k)
            forall (i=k+1:n, j=k:n+1) M(j,i)=M(j,i)-M(k,i)*M(j,k)
         enddo
         M(n+1,n)=M(n+1,n)/M(n,n)
         call solution(M,X,method)

         end subroutine gauss

         subroutine jordan(M,X)
         implicit none
         integer(mp) :: n,i,j,k
         real(mp) :: eps
         real(mp),dimension(1:) :: X(:)
         real(mp),dimension(1:size(X)+1,1:size(X)) :: M(:,:)
         character(10) :: method='jordan'
         eps=1.0_mp**(-mp)
         n = size(X, dim=1)         
         do k=1,n
            if (abs(M(k,k)) < eps) write(*,*) 'Ведущий элемент близок к нулю, номер строки:', k 
            forall (j=k:n+1) M(j,k)=M(j,k)/M(k,k)
            forall (i=1:k-1, j=k:n+1) M(j,i)=M(j,i)-M(k,i)*M(j,k)
            forall (i=k+1:n, j=k:n+1) M(j,i)=M(j,i)-M(k,i)*M(j,k)
         enddo
         M(n+1,n)=M(n+1,n)/M(n,n)
         call solution(M,X,method)

         end subroutine jordan

         subroutine lead(M,X)
         implicit none
         real(mp),dimension(1:)  ::  X(:)
         real(mp),dimension(1:size(X)+1,1:size(X)) :: M(:,:)         
         real(mp),allocatable  :: Just_array(:)
         integer(4),allocatable  :: Coordinate(:,:)
         integer(4) :: n,i,j,k
         character(10) :: method='lead'
         n = size(X, dim=1)
         allocate(Just_array(n+1))
         allocate(Coordinate(2,n))
         do k=1,n-1
            if ( abs( maxval(M(k:n,k:n)) ) > abs( minval(M(k:n,k:n)) ) ) then
              Coordinate(1:2,k)=(k-1)+maxloc(M(k:n,k:n))
            else
              Coordinate(1:2,k)=(k-1)+minloc(M(k:n,k:n))
            endif
              Just_array(1:n+1)=M(1:n+1,Coordinate(2,k))
              M(1:n+1,Coordinate(2,k))=M(1:n+1,k)
              M(1:n+1,k)=Just_array(1:n+1)
              Just_array(1:n)=M(Coordinate(1,k),1:n)
              M(Coordinate(1,k),1:n)=M(k,1:n)
              M(k,1:n)=Just_array(1:n)
              forall (j=k:n+1) M(j,k)=M(j,k)/M(k,k) 
              forall (i=k+1:n, j=k:n+1) M(j,i)=M(j,i)-M(k,i)*M(j,k)
         enddo
              M(n:n+1,n)=M(n:n+1,n)/M(n,n)
            call solution(M,X,method)
         do k=n-1,1,-1
              Just_array(k)=X(Coordinate(1,k))
              X(Coordinate(1,k))=X(k)
              X(k)=Just_array(k)
         enddo

         end subroutine lead

         subroutine solution(M,X,method)
         integer :: i, n
         real(mp),dimension(1:)  :: X(:)
         real(mp),dimension(1:size(X)+1,1:size(X)) :: M(:,:)         
         character(10) :: method
         n = size(X, dim=1)
         if (method /= 'jordan') then
            do i=n,1,-1
               X(i)=M(n+1,i)-dot_product(M(i+1:n,i),X(i+1:n))
            enddo
         else
         do i=1,n
               X(i)=M(n+1,i)
         enddo
         endif
         end subroutine solution

end module SLAU_solver