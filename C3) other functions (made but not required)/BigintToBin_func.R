
#Function to convert larger integers to binary (larger than R's intToBin function can perform)
#works only up to 2^66 which is big enough for our purposes but would need to be larger for networks eg. (66+ patches)

BigintToBin<-function(x, digits){
  x<-format(x, scientific=F)
  x<-as.character(x) #convert x to a character string
  x<-strsplit(x,"") #convert x to a vector of individual digits
  x<-lapply(x, FUN=as.numeric)
  x<-x[[1]]
  binary.x<-c() #create a vector to store each binary digit
  if (sum(x==0)==length(x)){ #if x is zero and therefor just "0"
    binary.x<-0 #in binary it's just 0
  }else{ #otherwise
    #while(sum(x[[1]]!=c("0"))!=0){ #until x becomes "0"
    while(sum(x==0)!=length(x)){
      if(odd(x[length(x)])){ #check if last digit of x is odd, if yes,
        binary.x<-c(binary.x, 1) #give a "1" as the next binary digit
      }else{ #otherwise,
        binary.x<-c(binary.x, 0) #give a "0" as the next binary digit
      }
      #then divide x by two, digit by digit left to right using the following algorithm:
      add<-0 #set the additive to be 0
      for(i in 1:length(x)){ #for each digit in x
        if (odd(x[i])){ #if the digit is odd
          next_add<-5 #set the next additive to be 5
        }else{ #otherwise,
          next_add<-0 #set the next additive to be 0
        }
        x[i]<-floor(x[i]/2)+add #replace that digit of x with 
        #the digit divided by two rounded to the next lowest integer and add the additive
        #and convert it back to a character
        add<-next_add #set the additive to be the next additive
        #if (x[[1]][i]=="0"){
        #  x[[1]]<-x[[1]][-i]
        #}
      }
      if(x[1]==0){ #remove 0 from the start of x if it has become a 0
        x<-x[-1]
      }
    }
  }
  if(digits!="NA"){ #if want the binary number given by a specific number of digits
    n.zeros<-digits-length(binary.x) #calculate how many zeros need to be added on
    if(n.zeros>=1){
      zeros<-rep(0, n.zeros)
      binary.x<-c(zeros, binary.x) #add the zero digits on
    }
  }
  return(rev(binary.x))
}

#TESTING:
bin<-BigintToBin((2^1), "NA") #WORKS! WOOO! takes only a couple secconds
length(bin)<(2^1)
bin<-BigintToBin((2^2), "NA") #WORKS! WOOO! takes only a couple secconds
length(bin)<(2^2)
bin<-BigintToBin((2^5), "NA") #WORKS! WOOO! takes only a couple secconds
length(bin)<(2^5)
bin<-BigintToBin((2^10), "NA") #WORKS! WOOO! takes only a couple secconds
length(bin)<(2^10)
bin<-BigintToBin((2^50), "NA") #WORKS! WOOO! takes only a couple secconds
length(bin)<(2^50)
#bin<-BigintToBin((2^100), "NA") #WORKS! WOOO! takes only a couple secconds
#length(bin)<(2^100)
#bin<-BigintToBin((2.8746218745621986428916*10^30), "NA") #WORKS! WOOO! takes only a couple secconds
#(2.8746218745621986428916*10^30)<(2^102)
#length(bin)<(2^102)
#bin<-BigintToBin((2.8746218745621986428916*10^30), (2^102)) #WORKS! WOOO! takes only a couple secconds


#rep(0, 2^50) #doesn't work
#rep(0, 2^30) #works
#rep(0, 2^35) #Starts but runs out of memory

#can quickly compute numbers in the nonillions into binary!!! :D
#BigintToBin((2^1000), "NA") #1000 takes too long for me to bother running
#BigintToBin(4, 4)
#x<-1576
#test1<-strsplit(intToBin(x),"")
#test1<-as.numeric(test1[[1]])
#test2<-BigintToBin(x,11)
#test1==test2

