#include <math.h>
#include <stdio.h>

typedef struct
{
    double x;
    double y;
    double z;
} Vector3D; // 三维向量

typedef struct
{
    Vector3D position;
    Vector3D velocity;
} Particlestate;

void get_electromagneticfield(Vector3D *E, Vector3D *B, double dipole_moment,
                              const Vector3D *position)
{
    double r = sqrt(position->x * position->x + position->y * position->y +
                    position->z * position->z);
    double pi = M_PI;
    double Vacuum_permeability = 4 * pi * 1e-7;

    E->x = 0;
    E->y = 0;
    E->z = 0;
    B->x = ((3 * Vacuum_permeability * dipole_moment) / (4 * pi * pow(r, 5))) *
           position->x * position->z;
    B->y = ((3 * Vacuum_permeability * dipole_moment) / (4 * pi * pow(r, 5))) *
           position->y * position->z;
    B->z = ((3 * Vacuum_permeability * dipole_moment) / (4 * pi * pow(r, 5))) *
           (3 * pow(position->z, 2) - r * r);
}

// 计算导数
void compute_derivative(const Particlestate *current, Particlestate *derivate,
                        double mass, double charge, double dipole_moment)
{
    Vector3D E, B;

    get_electromagneticfield(&E, &B, dipole_moment, &current->position);

    Vector3D v_cross_B;
    v_cross_B.x = current->velocity.y * B.z - current->velocity.z * B.y;
    v_cross_B.y = current->velocity.z * B.x - current->velocity.x * B.z;
    v_cross_B.z = current->velocity.x * B.y - current->velocity.y * B.x;

    double q_m = charge / mass;

    derivate->velocity.x = q_m * (E.x + v_cross_B.x);
    derivate->velocity.y = q_m * (E.y + v_cross_B.y);
    derivate->velocity.z = q_m * (E.z + v_cross_B.z);

    derivate->position = current->velocity;
}

void RK_method(Particlestate *state, double charge, double mass,
               double dipole_moment, double dt)
{
    Particlestate K1;
    Particlestate K2;
    Particlestate K3;
    Particlestate K4;
    Particlestate temp_state;

    // compute K1
    compute_derivative(state, &K1, mass, charge, dipole_moment);

    // compute K2
    temp_state = *state;
    temp_state.position.x += (dt / 2) * K1.position.x;
    temp_state.position.y += (dt / 2) * K1.position.y;
    temp_state.position.z += (dt / 2) * K1.position.z;
    temp_state.velocity.x += (dt / 2) * K1.velocity.x;
    temp_state.velocity.y += (dt / 2) * K1.velocity.y;
    temp_state.velocity.z += (dt / 2) * K1.velocity.z;
    compute_derivative(&temp_state, &K2, mass, charge, dipole_moment);

    // compute K3
    temp_state = *state;
    temp_state.position.x += (dt / 2) * K2.position.x;
    temp_state.position.y += (dt / 2) * K2.position.y;
    temp_state.position.z += (dt / 2) * K2.position.z;
    temp_state.velocity.x += (dt / 2) * K2.velocity.x;
    temp_state.velocity.y += (dt / 2) * K2.velocity.y;
    temp_state.velocity.z += (dt / 2) * K2.velocity.z;
    compute_derivative(&temp_state, &K3, mass, charge, dipole_moment);

    // compute K4
    temp_state = *state;
    temp_state.position.x += (dt)*K3.position.x;
    temp_state.position.y += (dt)*K3.position.y;
    temp_state.position.z += (dt)*K3.position.z;
    temp_state.velocity.x += (dt)*K3.velocity.x;
    temp_state.velocity.y += (dt)*K3.velocity.y;
    temp_state.velocity.z += (dt)*K3.velocity.z;

    compute_derivative(&temp_state, &K4, mass, charge, dipole_moment);

    state->position.x += (dt / 6) * (K1.position.x + 2 * K2.position.x +
                                     2 * K3.position.x + K4.position.x);
    state->position.y += (dt / 6) * (K1.position.y + 2 * K2.position.y +
                                     2 * K3.position.y + K4.position.y);
    state->position.z += (dt / 6) * (K1.position.z + 2 * K2.position.z +
                                     2 * K3.position.z + K4.position.z);
    state->velocity.x += (dt / 6) * (K1.velocity.x + 2 * K2.velocity.x +
                                     2 * K3.velocity.x + K4.velocity.x);
    state->velocity.y += (dt / 6) * (K1.velocity.y + 2 * K2.velocity.y +
                                     2 * K3.velocity.y + K4.velocity.y);
    state->velocity.z += (dt / 6) * (K1.velocity.z + 2 * K2.velocity.z +
                                     2 * K3.velocity.z + K4.velocity.z);
}

int main()
{
    double charge = 1.602e-19;
    double mass = 1.6726e-27;
    double dipole_moment = 7.94e22;
    double eV = 1.602e-19;
    double velocity = sqrt((2 * 1000000 * eV) / mass);
    Particlestate state = {
        {5 * (sqrt(2) / 2) * 6.371e6, 5 * (sqrt(2) / 2) * 6.371e6, 0},
        {velocity * cos(M_PI / 4), velocity * cos(M_PI / 3),
         velocity * sin(M_PI / 3)}};
    double dt = 5e-4;
    double T = 6e4;
    long long steps = T / dt;
    long long save_interval = 1e6;
    printf("总步数: %lld\n", steps);

    FILE *fp = fopen("track_RK.bin", "wb");
    if (fp == NULL)
    {
        printf("无法打开文件!\n");
        return 1;
    }

    long long save_interval_small = 100;   // 初始 100 步存一次
    long long save_interval_large = 50000; // 后期 50000 步存一次
    long long transition_step = 1e6;       // 在 1e6 步以后再降低存储密度
    double buffer[3];
    for (long long i = 0; i < steps; i++)
    {
        long long current_interval =
            (i < transition_step) ? save_interval_small : save_interval_large;

        if (i % current_interval == 0)
        {
            buffer[0] = state.position.x;
            buffer[1] = state.position.y;
            buffer[2] = state.position.z;
            fwrite(buffer, sizeof(double), 3, fp);
        }
        RK_method(&state, charge, mass, dipole_moment, dt);
    }

    fclose(fp);
    printf("计算完成，轨迹已保存到 track_RK.bin\n");

    return 0;
}
