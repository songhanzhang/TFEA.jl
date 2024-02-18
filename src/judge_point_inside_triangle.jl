function judge_point_inside_triangle(x1,y1,x2,y2,x3,y3,xp,yp)

    S1 = 1/2 * ( x3*yp + xp*y2 + x2*y3 - x3*y2 - x2*yp - xp*y3 )
    S2 = 1/2 * ( x1*yp + xp*y3 + x3*y1 - x1*y3 - x3*yp - xp*y1 )
    S3 = 1/2 * ( x2*yp + xp*y1 + x1*y2 - x2*y1 - x1*yp - xp*y2 )
    
    if S1 > 0 && S2 > 0 && S3 > 0
        answer = true
    else
        answer = false
    end
    
    return answer

end
