 
    module material_module
    use parameter_module
    
    implicit none
  
    
    
    type :: material ! global material array is on this type
        
        private
        
        character(len=matnamelength) :: matname=''  ! user-given name from input
        character(len=mattypelength) :: mattype=''  ! defined material types: isotropic, lamina, interlayer, etc.
        integer :: typekey=0                         ! index in the respective array of the associated material type
        
    end type
    
    
    

    interface empty
        module procedure empty_material
    end interface
    
    interface update
        module procedure update_material
    end interface
    
    interface extract
        module procedure extract_material
    end interface
    
    private :: empty_material, update_material, extract_material
    
    
    
    contains
    
    
    
    
    ! this subroutine is used in the preprocessing to format the 
    ! global material library lib_mat
    subroutine empty_material(mat)
        type(material),intent(out):: mat
        
        mat%matname=''
        mat%mattype=''
        mat%typekey=0
    
    end subroutine empty_material
    
    
    ! this subroutine is used in the preprocessing to fill in the
    ! material information in the global material library lib_mat
    subroutine update_material(mat,matname,mattype,typekey)
        type(material),intent(inout):: mat
        character(len=*),optional,intent(in) :: matname 
        character(len=*),optional,intent(in) :: mattype 
        integer,optional,intent(in) :: typekey
        
        if(present(matname)) mat%matname=matname
        if(present(mattype)) mat%mattype=mattype
        if(present(typekey)) mat%typekey=typekey
    
    end subroutine update_material
    
    
    ! this subroutine is used anywhere to extract material information
    subroutine extract_material(mat,matname,mattype,typekey)
        type(material),intent(in):: mat
        character(len=matnamelength),optional,intent(out) :: matname 
        character(len=mattypelength),optional,intent(out) :: mattype 
        integer,optional,intent(out) :: typekey
        
        if(present(matname)) matname=mat%matname
        if(present(mattype)) mattype=mat%mattype
        if(present(typekey)) typekey=mat%typekey
    
    end subroutine extract_material
    
    
    end module material_module